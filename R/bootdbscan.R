####################
### bootdbscan() ###
####################
#
#' @title DBSCAN for big data using a bootstrap approach
#'
#' @description
#' The function samples observations and performs DBSCAN clustering (based on [dbscan::dbscan()]) on subsamples. Clusters are denoted binary (1/0) as cluster/noise. Results are aggregated
#' by observations and majority voting is used to flag an observation as cluster or noise.
#'
#' @details
#' The function runs in a loop DBSCAN algorithm (based on [dbscan::dbscan()]) using subsamples smaller than a full sample. Clusters are binary coded (cluster/noise) and aggregated by observation ID.
#' For each row (observation), a majority vote indicates whether a point is classified as cluster or noise. As a rule of thumb, `sample_size*times=150%` to `200%` of the original
#' number of observations. Points that were never selected for subsampling are classified as cluster/noise using the k nearest neighbours algorithm - they are given
#' the same flag as their first nearest neighbour.
#'
#' @name bootdbscan
#' @param data_sf Object in `sf` or the `data.frame` class - in the case of a `data.frame` object, the first and second columns must contain X and Y coordinates.
#' @param sample_size The size of the subset used in a single clustering run.
#' @param times The number of times clustering is run on a random sub-sample.
#' @param eps [dbscan::dbscan()] parameter to define a radius around each point to count neighbouring points.
#' @param minPts The minimum number of points (in a radius eps) to classify the surroundings as dense, [dbscan::dbscan()] parameter.
#' @param plot Logical; indicates whether the function should generate a plot, default `plot=TRUE`.
#'
#' @return `bootdbscan()` returns an object of class `bootdbscan` with the following components:
#' \item{cluster }{An integer vector of cluster assignments. 0 indicates noise points. 1 indicates high density cluster assignments.}
#' \item{sample_size }{ value of the `sample_size` parameter.}
#' \item{times }{ value of the `times` parameter.}
#' \item{eps }{ value of the `eps` parameter.}
#' \item{minPts }{ value of the `minPts` parameter.}
#'
#' @examples
#' bdb<-bootdbscan(firms_sf, sample_size=1000, times=20, eps=0.1, minPts=5)
#'
#' @export
bootdbscan<-function(data_sf, sample_size, times, eps, minPts, plot=TRUE){
  if((inherits(data_sf,"sf") && st_geometry_type(data_sf,FALSE)=="POINT")) {
    crds<-as.data.frame(st_coordinates(data_sf))
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class sf.\n", sep = "")
  }  else if(inherits(data_sf,"data.frame",TRUE)==1){
    crds<-data_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class data.frame.\n", sep = "")
  }  else {
    stop("The class of data_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # zbadać sample_size i ustawić ew. na wielkość zbioru danych
  sample_size <- as.integer(sample_size)
  if (sample_size > nrow(crds) || sample_size<1) {
    sample_size<-nrow(crds)
    cat("Wrong sample size. Sample_size set to:",sample_size,"\n",sep="")
  }

  # zbadać times, czy nie ujemne - może inny warunek niż <=0
  times<-round(times,0)
  if (times<=0) {
    times<-16
    cat("Wrong number of clusters. Number of cluster set to the proposed value: 16.","\n",sep="")
  }

  selector<-matrix(0, nrow=sample_size, ncol=times)
  result<-matrix(0, nrow=sample_size, ncol=times)

  for(i in 1:times) selector[,i]<-sample(1:nrow(crds), sample_size, replace=FALSE) # sample points

  # poniższe można uprościć (trochę) - zmniejszyć liczbę pośrednich obiektów
  for(i in 1:times){
    vec<-selector[,i]
    sub<-as.matrix(crds[vec,])
    dbs<-dbscan(sub, eps=eps, minPts=minPts) # dbscan in subgroups
    result[,i]<-dbs$cluster}

  DBS<-data.frame(ID=1:nrow(crds), new_cluster=rep(NA, times=nrow(crds)))

  selector_v<-as.vector(selector)
  result_v<-as.vector(result)
  result_v<-ifelse(result_v>0, 1, 0)
  x<-cbind(selector_v, result_v)

  x.grouped <- summarise(group_by(as.data.frame(x), selector_v), mean = mean(result_v))

  DBS$new_cluster[x.grouped$selector_v]<-round(x.grouped$mean,0) # for the points that were sampled at least once,
  # assign 0 if at average they were outliers, and 1
  # if at average they were in a cluster.

  no_cluster<-which(is.na(DBS$new_cluster))
  nneigh = 1

  while(length(no_cluster)>0){  # as long as there are points with no assignment...
    crds.knn<-kNN(crds, nneigh, sort=FALSE) # ...find their next nearest neighbour...
    DBS$new_cluster[no_cluster]<-DBS$new_cluster[crds.knn$id[no_cluster,nneigh]] #...and assign the assignment of the neighbour
    nneigh = nneigh+1 # check by the further neighbour
    no_cluster<-which(is.na(DBS$new_cluster))}

  # wykres do sprawdzenia czy o to chodziło, tytuł do sprawdzenia
  if (plot){
    plot(crds,pch=".", main="Observations with core points (red)")
    points(crds[DBS$new_cluster==1,],pch=".",col="red")
  }

  structure(
    list(
      cluster = DBS$new_cluster,
      sample_size = sample_size,
      times = times,
      eps = eps,
      minPts = minPts
    ),
    class = "bootdbscan"
  )

}

#' @export
# Define the output format for print(bootdbscan)
print.bootdbscan <- function(x, ...) {
  cat("BOOTDBSCAN clustering for", length(x$cluster), "objects.\n")
  cat("Parameters: sample_size = ", x$sample_size,", times = ",x$times, ", eps = ", x$eps, ", minPts = ", x$minPts, "\n", sep="")
  cat("The clustering contains ", length(x$cluster[x$cluster == 1]), "cluster points and ", length(x$cluster[x$cluster == 0]), " noise points.\n")
  cat("Available fields: cluster, sample_size, times, eps, minPts\n")
}

#' @rdname bootdbscan
#' @export
BOOTDBSCAN <- bootdbscan
