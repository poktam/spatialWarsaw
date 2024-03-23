####################
### ClustConti() ###
####################
#
#' @title Division of geolocalised points into spatially continuous clusters using DBSCAN
#'
#' @description
#' This is a wrapper function for the DBSCAN algorithm (based on [dbscan::dbscan()]). For a given number of n clusters, it generates the division into n+1 high-density groups
#' (n clusters + noise group) that guarantees the noise ratio specified by the user.
#'
#' @details
#' The function performs modified DBSCAN clustering (based on [dbscan::dbscan()]), with user-specified restrictions on the number of clusters and noise ratio. It produces spatially
#' continuous clusters of geolocated points. The resulting clusters may have different internal densities and sizes. This method assumes that the area
#' and its spatial structure are known (as visible agglomeration points), which allows an appropriate number of clusters to be set. For example,
#' for urbanisation data, each cluster will cover the whole city with suburbs.
#'
#' The method adjusts the basic parameters of DBSCAN: minPts and eps, in iterations, to obtain user-defined parameters of the final partitioning:
#' ratio of noise and number of clusters.
#'
#' The output of this algorithm is identical to typical DBSCAN, but the input parameters are different (number of clusters and percentage of noise,
#' instead of the number of points within the radius minPts and radius eps) - this is to get control over the relative point density and the number of clusters.
#'
#' @name ClustConti
#' @param data_sf Object in `sf` or the `data.frame` class - in the case of a `data.frame` object, the first and second columns must contain X and Y coordinates.
#' @param clusters Number of density clusters to generate; function automatically adds an extra cluster for noise data (low density cluster).
#' @param noise Proportion of noise data that should result from DBSCAN clustering, can be between 0 (no noise) and 1 (all noise).
#' @param r_p DBSCAN optimisation parameter for radius, default value 0.001.
#' @param eps_r DBSCAN parameter - the initial size of the radius in coordinate units, default value 1e-15.
#' @param eps_np DSBCAN optimisation parameter for radius, default value 0.01.
#' @param minPts0 DBSCAN parameter - expected number of points within a radius to be considered as being dense, default value 5.
#'
#' @return `ClustConti()` returns a list object:
#' \item{type}{Simply the name of the wrapper function used: "ClustConti".}
#' \item{cluster}{Vector of cluster assignments, 1:n are high density clusters, 0 is noise (low density cluster of remaining points).}
#' \item{parameters}{A list of 6 parameters used in clustering - two user-defined (`clusters` and `noise`) and four that can also be set
#' by default (`r_p`, `eps_r`, `eps_np`, `minPts0`).}
#'
#' @seealso [ClustDisjoint()], [ssr()]
#'
#' @examples #To be done!!!
#'
#' @export
ClustConti<-function(data_sf, clusters, noise, r_p=0.001, eps_r=10e-16, eps_np=10e-3, minPts0=5){
  params <- as.list(environment())[-1]

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

  # zbadać clusters, czy nie ujemne - może inny warunek niż <=0
  clusters<-round(clusters,0)
  if (clusters<=0) {
    clusters<-4
    cat("Incorrect clusters parameter. Set to the proposed value: 4.","\n",sep="")
  }

  # zbadać minPts0, czy nie ujemne - może inny warunek niż <=0
  minPts0<-round(minPts0,0)
  if (minPts0<=0) {
    minPts0<-5
    cat("Incorrect minPts0 parameter. Set to the proposed value: 5.","\n",sep="")
  }

  minPts<-minPts0
  sink(nullfile()) #ew. spróbować z capture.output() i invisible()
  output<-ClustDisjoint(data_sf=crds, noise=noise, r_p=r_p, eps_r=eps_r, eps_np=eps_np, minPts=minPts)$cluster # run ClustDisjoint with initial parameters
  sink()
  minPts_prev<-minPts

  if(length(unique(output))>clusters+1) {
    minPts=2*minPts
  } else {
    minPts=round(minPts/2,0)} # adjust minPts
  output_prev<-output

  while(length(unique(output))!=clusters+1 & minPts_prev!=minPts){ # repeat below until the desired number of clusters is obtained
    minPts=minPts+1

    sink(nullfile())
    output<-ClustDisjoint(data_sf=crds, noise=noise, r_p=r_p, eps_r=eps_r, eps_np=eps_np, minPts=minPts)$cluster # run ClustDisjoint
    sink()

    if(length(unique(output))>clusters+1 & length(unique(output_prev))>clusters+1) {
      minPts_prev=minPts
      minPts=2*minPts } else if(length(unique(output))<clusters+1 & length(unique(output_prev))<clusters+1){
        minPts_prev=minPts
        minPts=round(minPts/2,0)} else {
          rr=minPts
          minPts=round((minPts+minPts_prev)/2,0)
          minPts_prev=rr  # adjust minPts
        }
    output_prev<-output

  }
  cat("\nFinal value: minPts = ", minPts,"\n",sep="")
  if(length(unique(output))!=clusters+1) warning(paste0("Failed to get ", clusters, " clusters. Try changing minPts0 value."))

  # zastanowić się, czy klasa testing nie jest do usunięcia?
  structure(
    list(
      type = "ClustConti",
      cluster = output,
      parameters = params
    ),
    class = "testing"
  )

}


#######################
### ClustDisjoint() ###
#######################
#
#' @title Division of geolocalised points into spatially disjoint clusters using DBSCAN
#'
#' @description
#' This is a wrapper function for the DBSCAN algorithm (based on [dbscan::dbscan()]). For a given number of n clusters (given as n-1 interval noise thresholds),
#' it generates the division into n spatially disjoint density groups that guarantee the similar density of points within each cluster.
#' It also works in the bootstrap version.
#'
#' @details
#' The function performs modified DBSCAN clustering (based on [dbscan::dbscan()]), which generates `n` spatially disjoint clusters of similar density from geolocated
#' points - the result is like concentric density rings around cities. The resulting clusters have similar internal densities and similar sizes.
#' The clustering depends on two user-specified parameters: `minPts` - the required number of points within the radius (argument typically used
#' in DBSCAN clustering), and `noise` - a vector of `n-1` density thresholds sorted in descending order: `noise=(𝑠1, 𝑠2, ..., 𝑠(n-1)`) such
#' that `s1>𝑠2>...>𝑠(n-1)`, specifying the desired density levels in `n` groups, measured as a percentage of points classified as noise
#' in each iteration of the algorithm. The number of clusters is implicitly specified - for `n` clusters, the user specifies `n-1` thresholds.
#'
#' The algorithm iteratively labels points as clustered or noise - in the first step it finds the most dense locations so that `(1-s1)%` of points
#' are classified as core, `minPts` is fixed and radius eps is adjusted. In the second step, it considers all points again and repeats the clustering
#' into core or noise groups, with the same rules and with a lower noise threshold `s2` - it classifies points into group 1 (from step 1),
#' group 2 (from step 2) or as noise. Clustering stops when there are `n-1` groups of core points and the `Nth` group of noise (including `𝑠(n-1)%` of points).
#' The solution is highly dependent on the number of `minPts`, the rule of thumb to set `minPts=1/100` of the points works well.
#'
#' @name ClustDisjoint
#' @param data_sf Object in `sf` or the `data.frame` class - in the case of a `data.frame` object, the first and second columns must contain X and Y coordinates.
#' @param noise A vector of n-1 elements for n clusters, sorted in descending order; proportion of noise data in each interval, e.g. c(0.8, 0.6, 0.4, 0.2)
#' for n=5 clusters of similar internal density.
#' @param r_p DBSCAN optimisation parameter for radius, default value 0.001.
#' @param eps_r DBSCAN parameter - the initial size of the radius in coordinate units, default value 1e-15.
#' @param eps_np DSBCAN optimisation parameter for radius, default value 0.01.
#' @param minPts DBSCAN parameter - expected number of points within a radius to be considered as being dense, devault value 5
#' @param bootstrap Logical; `FALSE` for standard approach, `TRUE` for bootstrap approach (more in the help for [bootdbscan()]),
#' suitable for big data sets.
#' @param sample_size The number of observations used in each bootstrap iteration should be less than the total number of observations,
#' used only if `bootstrap=TRUE`.
#' @param times The number of iterations in the bootstrap, we recommend that `sample_size*times = 150%` or `200%` of the original number of observations,
#' used only when `bootstrap=TRUE`.
#'
#' @return `ClustDisjoint()` returns a list object:
#' \item{type}{Simply the name of the wrapper function used: "ClustDisjoint".}
#' \item{cluster}{Vector of cluster assignments, 1:n-1 are high density clusters, n is noise (low density cluster of remaining points).}
#' \item{parameters}{A list of 8 parameters used in clustering - `noise`, `r_p`, `eps_r`, `eps_np`, `minPts`, `bootstrap`,`sample_size`, `times`.}
#'
#' @seealso [ClustConti()], [ssr()]
#'
#' @examples #To be done!!!
#'
#' @export
ClustDisjoint<-function(data_sf, noise=c(0.8, 0.6, 0.4, 0.2), r_p=0.001, eps_r=10e-16, eps_np=10e-3, minPts=5,
                bootstrap=FALSE, sample_size=NULL, times=NULL){
  params <- (as.list(environment()))[-1] # tu można poprawić, żeby w zależności od bootstrap dawało odpowiednie

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

  # zbadać minPts, czy nie ujemne - może inny warunek niż <=0
  minPts<-round(minPts,0)
  if (minPts<=0) {
    minPts<-5
    cat("Incorrect minPts parameter. Set to the proposed value: 5.","\n",sep="")
  }

  if(bootstrap){

    if (is.null(sample_size) || is.null(times)) {
      stop("For bootstrap=TRUE, both sample_size and times parameters must be set.")
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
      cat("Incorrect times value. Times set to the proposed value: 16.","\n",sep="")
    }
  }

  noise<-sort(noise, decreasing=TRUE)
  output<-rep(length(noise)+1,nrow(crds))

  r<-r_p

  # można będzie spróbować zmniejszyć kod w dużym if (poniżej), bo się powtarza
  if(bootstrap){
    for(i in 1:length(noise)) {
      sink(nullfile())
      result_temp<-bootdbscan(data_sf=crds, sample_size=sample_size, times=times, eps=r, minPts=minPts, plot=FALSE)$cluster # run bootdbscan with initial parameters
      sink()
      noisePercent<-length(result_temp[result_temp==0])/nrow(crds)
      r_prev=r

      if(noisePercent>noise[i]) r=2*r else r=r/2 # adjust the radius
      nP_prev<-noisePercent

      while((abs(r-r_prev)>eps_r) & (abs(noisePercent-noise[i])>eps_np)){ # repeat below until the desired percentage of noise is obtained
        sink(nullfile())
        result_temp<-bootdbscan(data_sf=crds, sample_size=sample_size, times=times, eps=r, minPts=minPts, plot=FALSE)$cluster # run bootdbscan
        sink()
        noisePercent<-length(result_temp[result_temp==0])/nrow(crds)

        if(noisePercent>noise[i] & nP_prev>noise[i]) { # adjust theradius
          r_prev=r
          r=2*r } else if(noisePercent<noise[i] & nP_prev<noise[i]){
            r_prev=r
            r=r/2} else {
              rr=r
              r=(r+r_prev)/2
              r_prev=rr}

        nP_prev<-noisePercent

      }

      output<-ifelse(output==length(noise)+1 & result_temp>0, i, output) # final result
    }

  }else {

    for(i in 1:length(noise)) {
      result_temp<-dbscan(crds, r, minPts=minPts) # run dbscan
      noisePercent<-length(result_temp$cluster[result_temp$cluster==0])/nrow(crds)
      r_prev=r

      if(noisePercent>noise[i]) r=2*r else r=r/2 # adjust the radius
      nP_prev<-noisePercent

      while(abs(r-r_prev)>eps_r & abs(noisePercent-noise[i])>eps_np){ # repeat the following until the desired percentage of noise is obtained
        result_temp<-dbscan(crds, r, minPts=minPts) # run dbscan
        noisePercent<-length(result_temp$cluster[result_temp$cluster==0])/nrow(crds)

        if(noisePercent>noise[i] & nP_prev>noise[i]) { # adjust the radius
          r_prev=r
          r=2*r } else if(noisePercent<noise[i] & nP_prev<noise[i]){
            r_prev=r
            r=r/2} else {
              rr=r
              r=(r+r_prev)/2
              r_prev=rr}

        nP_prev<-noisePercent
      }

      output<-ifelse(output==length(noise)+1 & result_temp$cluster>0, i, output)} # adjust the radius

    if(length(noise)==1) output<-result_temp$cluster}

  # zastanowić się, czy klasa testing nie jest do usunięcia?
  structure(
    list(
      type = "ClustDisjoint",
      cluster = output,
      parameters = params
    ),
    class = "testing"
  )

}

# Define the output format for print(testing)
#' @export
print.testing <- function(x) {
  cat("Type", x$type,"group assignments for", length(x$cluster), "objects.\n")
  cat("Summary of assignments:\n")
  print(summary(factor(x$cluster)))
  cat("\nAvailable fields: type, cluster, parameters\n")
}
