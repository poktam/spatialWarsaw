#############
### QDC() ###
#############
#
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name QDC
#' @aliases QDC qdc
#' @param data_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords)
#' @param sample_size Sample size, must be less than or equal to the number of points in the dataset (`data_sf` parameter). If `sample_size` is larger, it is automatically set to the number of points in the dataset.
#' @param nclust Number of clusters to find. Must be >=2. For nclust=3 additional cluster labels will be included in the output object.
#' @param k Number of nearest neighbors used for the distance calculation. Default value set to 10.
#' @param eps do opisu
#' @examples #To be done!!!
#'
#' @return `QDC()` returns ... to be done.
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
  ggplot(data.d, aes(x=knndist1.scaled, y=frnn1.scaled, color=km.set1)) +  geom_point()+ theme(legend.position="none") # xy plot of spat.var

  # location of clusters in space
  ggplot(data.d, aes(x=X, y=Y, color=km.set1)) +  geom_point()+ theme(legend.position="none")

  # LEPSZE NAZWY KOLUMN W OUTPUCIE
  # czy potrzeba zwracać listę parametrów i obiekt output czy tylko obiekt output wystarczy?
  # poprawić, żeby może na ekranie było jakieś summary/table a nie lista
  return(data.d)
}



