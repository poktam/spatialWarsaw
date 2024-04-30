#############
### QDC() ###
#############
#
#
#' @title Quick Density Clustering (QDC): identifying clusters of different densities for geolocated points
#'
#' @description
#' Quick Density Clustering (QDC) uses two spatial variables - the total distance to k nearest neighbours
#' and the number of neighbours within a fixed radius (`eps`). Both spatial variables are normalised
#' and clustered using the K-means algorithm. Clusters separate points with different local densities.
#' The function allows sampling to speed up the calculations.
#'
#' @details
#' This algorithm is an alternative to DBSCAN. It has no problem with parameterisation, especially in the case of sampling -
#' it gives a robust solution. It uses K-means, but unlike most algorithms - it does not cluster spatial coordinates,
#' but pre-prepared normalised spatial variables: total distance to k nearest neighbours and number of neighbours
#' within a fixed radius (`eps`). Its logic is simple: densely located points have a low total distance to `knn` points,
#' while the number of points within a given radius is high. Sparsely located points, on the other hand, have a high total distance
#' to `knn` points, while the number of nearest neighbours within a given radius is low. Both variables are self-scaling
#' after normalisation and easily show what is high or low. K-means clustering naturally separates groups with similar values
#' of both spatial variables. There is a non-linear relationship between the two spatial variables, visible in the output Figure 1.
#' K-means clustering easily finds thresholds of both variables that separate clusters.
#'
#' The algorithm is robust to sampling. Limiting the size of the data set speeds up the computation. The parameter `sample_size` allows to find
#' the best relation between quality and computation time.
#'
#' @name QDC
#' @param data_sf Object in `sf` or the `data.frame` class - in the case of a `data.frame` object, the first and second columns must contain X and Y coordinates.
#' @param sample_size Sample size, must be less than or equal to the number of points in the dataset (`data_sf` parameter).
#' If `sample_size` is larger, it is automatically set to the number of points in the dataset.
#' @param nclust Number of clusters to find. Must be >=2. For nclust=3 additional cluster labels are included in the output object.
#' @param k Number of nearest neighbours used to calculate the total distance. The default value is k=10.
#' @param eps Radius for counting neighbours, fixed value for all observations. The default value is eps=0.05.
#' This is approximately 5 km for points in the WGS84 projection.
#'
#' @return `QDC()` returns two visualisations:
#' * Scatter plot of both spatial variables coloured as clusters
#' * Location of spatial points coloured according to cluster membership.
#'
#' `QDC()` returns also a data.frame object containing the following columns:
#' \item{X}{original X coordinates}
#' \item{Y}{original Y coordinates}
#' \item{knndist1}{total distance to k nearest neighbours}
#' \item{frnn1}{number of nearest neighbours with fixed radius}
#' \item{knndist1.scaled}{spatial variable knndist1 after normalisation}
#' \item{frnn1.scaled}{spatial variable frnn1 after normalisation}
#' \item{km.set1}{cluster ID of the observation, from the K-means algorithm}
#' \item{outcome.set1}{cluster labels, only if k=3}
#'
#' @references
#' Kopczewska K., (under review), QDC: Quick Density Clustering of Geo-located Data
#'
#' @examples #To be done!!!
#'
#' @export
QDC<-function(data_sf, sample_size, nclust=3, k=10, eps=0.05){
  if((inherits(data_sf,"sf") && st_geometry_type(data_sf,FALSE)=="POINT")) crds<-sf::st_coordinates(data_sf)
  else if(inherits(data_sf,"data.frame",TRUE)==1){
    crds<-data_sf[,c(1,2)]
    colnames(crds)<-c("X","Y")
  }
  else {
    stop("The class of data_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # zbadać sample_size i ustawić ew. na wielkość zbioru danych
  sample_size <- as.integer(sample_size)
  if (sample_size > nrow(crds) || sample_size<1) {
    sample_size<-nrow(crds)
    cat("Wrong sample size. Sample_size set to:",sample_size,"\n",sep="")
  }

  # zbadać n, czy nie ujemne
  nclust<-round(nclust,0)
  if (nclust<=2) {
    nclust<-3
    cat("Wrong number of clusters. Number of cluster set to the default value: 3.","\n",sep="")
  }

  # zbadać k, czy nie ujemne
  k<-round(k,0)
  if (k<1) {
    k<-10
    cat("Wrong number of nearest neighbors. Number of nearest neighbors set to the default value: 10.","\n",sep="")
  }

  # część wspólna niezależna od nclust
  data.d<-as.data.frame(crds[sample(nrow(crds), sample_size, replace=FALSE), ])

  # spatial variable – sum of distances to knn
  knn.dist<-kNNdist(as.matrix(data.d), k, all =TRUE)
  data.d$knndist1<-apply(knn.dist,1,sum)

  # spatial variable – number of points in fixed radius
  agg.radius<-frNN(as.matrix(data.d), eps=0.05)
  data.d$frnn1<-unlist(lapply(agg.radius$id, length))

  # normalization of spatial variables
  data.d$knndist1.scaled<-scale(data.d$knndist1)
  data.d$frnn1.scaled<-scale(data.d$frnn1)

  # k-means clustering of two normalized spatial variables
  km.set1<-kmeans(data.d[ ,5:6], nclust)
  data.d$km.set1<-as.factor(km.set1$cluster)

  # etykiety dla nclust=3, może jakoś to uprościć (thresholds)
  if (nclust==3) {
    # thresholds
    t1<-max(min(data.d$knndist1.scaled[data.d$km.set1==1]),
            min(data.d$knndist1.scaled[data.d$km.set1==2]),
            min(data.d$knndist1.scaled[data.d$km.set1==3])) # threshold
    t1	# when knndist>t1 – it is low-density cluster

    t2<-max(min(data.d$frnn1.scaled[data.d$km.set1==1]),
            min(data.d$frnn1.scaled[data.d$km.set1==2]),
            min(data.d$frnn1.scaled[data.d$km.set1==3]))  # threshold
    t2	# when frnn(agg)>t2 – it is high-density cluster

    # classification of points to clusters
    # ten krok jest tylko dla 3 klastrów, bo jest przewidywalne
    data.d$outcome.set1<-ifelse(data.d$knndist1.scaled>t1, "low-density", ifelse(data.d$frnn1.scaled>t2,"high-density", "mid-density"))
  }


  # visualization of clusters in 2D of spatial variables
  p1<-ggplot(data.d, aes_string(x="knndist1.scaled", y="frnn1.scaled", color="km.set1")) +  geom_point()+ theme(legend.position="none") # xy plot of spat.var

  # location of clusters in space
  p2<-ggplot(data.d, aes_string(x="X", y="Y", color="km.set1")) +  geom_point()+ theme(legend.position="none")

  grid.arrange(p1,p2, ncol = 2, nrow = 1)

#  print(p1)
#  print(p2)

  # LEPSZE NAZWY KOLUMN W OUTPUCIE
  # czy potrzeba zwracać listę parametrów i obiekt output czy tylko obiekt output wystarczy?
  # poprawić, żeby może na ekranie było jakieś summary/table a nie lista
  return(data.d)
}

#' @rdname QDC
#' @export
qdc <- QDC


