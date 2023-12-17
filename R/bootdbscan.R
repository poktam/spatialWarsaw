####################
### bootdbscan() ###
####################
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name bootdbscan
#' @aliases bootdbscan BOOTDBSCAN
#' @param crds do opisu
#' @param sample_size do opisu
#' @param times do opisu
#' @param eps do opisu
#' @param minPts do opisu
#' @examples #To be done!!!
#'
#' @return `bootdbscan()` returns an object of class `bootdbscan` with the following components:
#'
#' \item{cluster }{A integer vector with cluster assignments. Zero indicates noise points.}
#' \item{sample_size }{ value of the `sample_size` parameter.}
#' \item{times }{ value of the `times` parameter.}
#' \item{eps }{ value of the `eps` parameter.}
#' \item{minPts }{ value of the `minPts` parameter.}
#'
#' @export
bootdbscan<-function(crds, sample_size, times, eps, minPts){

  selector<-matrix(0, nrow=sample_size, ncol=times)
  result<-matrix(0, nrow=sample_size, ncol=times)

  for(i in 1:times) selector[,i]<-sample(1:dim(crds)[1], sample_size, replace=FALSE) # sample points

  for(i in 1:times){
    vec<-selector[,i]
    sub<-as.matrix(crds[vec,])
    dbs<-dbscan::dbscan(sub, eps=eps, minPts=minPts) # dbscan in subgroups
    result[,i]<-dbs$cluster}

  DBS<-data.frame(ID=1:dim(crds)[1], new_cluster=rep(NA, times=dim(crds)[1]))

  selector_v<-as.vector(selector)
  result_v<-as.vector(result)
  result_v<-ifelse(result_v>0, 1, 0)
  x<-cbind(selector_v, result_v)

  x.grouped <- dplyr::summarise(dplyr::group_by(as.data.frame(x), selector_v), mean = mean(result_v))

  DBS$new_cluster[x.grouped$selector_v]<-round(x.grouped$mean,0) # for the points that were sampled at least once,
  # assign 0 if at average they were outliers, and 1
  # if at average they were in a cluster.

  no_cluster<-which(is.na(DBS$new_cluster))
  nneigh = 1

  while(length(no_cluster)>0){  # as long as there are points with no assignment...
    crds.knn<-dbscan::kNN(crds, nneigh, sort=FALSE) # ...find their next nearest neighbour...
    DBS$new_cluster[no_cluster]<-DBS$new_cluster[crds.knn$id[no_cluster,nneigh]] #...and assign the assignment of the neighbour
    nneigh = nneigh+1 # check by the further neighbour
    no_cluster<-which(is.na(DBS$new_cluster))}

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
#' @export
BOOTDBSCAN <- bootdbscan
