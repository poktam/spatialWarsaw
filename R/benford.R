#########################
### SpatBenfordTest() ###
#########################
#
# SPRAWDZIĆ CZY NIE OPRZEĆ NA PAKIECIE benford (nowszy) zamiast benford.analysis
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name SpatBenfordTest
#' @aliases SpatBenfordTest spatbenfordtest
#' @param data_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords)
#' @param var_name Name of column with additional variable. If empty, a 2d distribution is tested, if given, a 3d distribution.
#' @param sample_size Sample size, must be less than or equal to the number of points in the dataset (`data _sf` parameter). If `sample_size` is larger, it is automatically set to the number of points in the dataset. We suggest that a value greater than 800 is not used for reasons of computational efficiency. (SPRAWDZIĆ)
#' @examples #To be done!!!
#'
#' @return `SpatBenfordTest()` returns ... to be done.
#'
#' @export
SpatBenfordTest<-function(data_sf, sample_size, var_name=NULL){
  # koordynaty w zależności od typu danych
  #DODAĆ SPRAWDZENIE CZY SF W FORMIE POINTS!!!!
  if(inherits(data_sf,"sf")) crds<-sf::st_coordinates(data_sf)
  else if(inherits(data_sf,"data.frame",TRUE)==1){
    crds<-data_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    }
  else {
      stop("The class of data_sf must be 'sf' or 'data.frame' only.")
  }

  # zbadać sample_size i ustawić ew. na wielkość zbioru danych
  sample_size <- as.integer(sample_size)
  if (sample_size > nrow(crds) || sample_size<1) {
    sample_size<-nrow(crds)
    cat("Wrong sample size. Sample_size set to:",sample_size,"\n",sep="")
  }

  # zbadać var_name - jeśli nie ma lub nie istnieje w danych, to pominąć, jeśli istnieje to wyciągnąć kolumnę
  if (length(var_name)>1) {
    var_col<-var_name[1]
    cat("Parameter var_name longer than 1. The first element has been selected: ",var_name,"\n",sep="")
  }

  m <- match(gsub(" ", ".", var_name), colnames(data_sf))

  if (is.na(m) || is.null(var_name)) {
    cat("Unknown variable name or variable name not specified. 2D spatial distribution will be tested.","\n",sep="")
  } else {
    crds<-cbind(crds,data_sf[,m])
    colnames(crds)[3]<-colnames(data_sf)[m]
    cat("3D spatial distribution will be tested.","\n",sep="")
  }

  # sample do testowania
  data.d<-crds[sample(nrow(crds), sample_size, replace=FALSE), ]

  my.dist<-dist(data.d)
  benford_result<-benford(as.vector(my.dist))
  plot(benford_result)
  return (benford_result)
}




###############################
### SpatBenfordPattern() ###
###############################
# R code for SpatBenfordPattern()
# generating spatial Benford-like 2D distribution
################################
#
# SPRAWDZIĆ CZY NIE OPRZEĆ NA PAKIECIE benford (nowszy) zamiast benford.analysis
# DODAĆ TESTOWANIE WYNIKU ZA POMOCĄ SPATBENFORDPATTERN()
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name SpatBenfordPattern
#' @aliases SpatBenfordPattern spatbenfordpattern
#' @param region_sf do opisu (obiekt sf ale jako region)
#' @param sample_size Number of points to be generated based on the region. Default value 5000. We suggest that a value greater than ..... is not used for reasons of computational efficiency. (SPRAWDZIĆ, UWAGA! st_sample czasem generuje mniej punktów)
#' @references #To be done - ważne tutaj bo sztywne parametry w kodzie
#' @examples #To be done!!!
#'
#' @return `SpatBenfordPattern()` returns ... to be done.
#'
#' @export
SpatBenfordPattern<-function(region_sf, sample_size=5000){
  #DODAĆ SPRAWDZENIE CZY SF W FORMIE REGIONU!!!!
  if(!inherits(region_sf,"sf")) {
    stop("The class of region_sf must be 'sf' only.")
  }

  ss_share<-c(0.810, 0.120, 0.070) # % of obs. clustered / regular /random
  bbox<-st_bbox(region_sf)

  # clustered skewed distribution
  x<-rnorm(sample_size*ss_share[1], bbox[1]+0.25*(bbox[3]-bbox[1]), 0.08*(bbox[3]-bbox[1]))
  y<-rnorm(sample_size*ss_share[1], bbox[2]+0.75*(bbox[4]-bbox[2]), 0.08*(bbox[4]-bbox[2]))
  bpp1<-data.frame(X=x, Y=y)

  bpp_sf2<-st_sample(region_sf, round(sample_size*ss_share[2],0), type="regular")
  bpp2<-st_coordinates(bpp_sf2)

  bpp_sf3<-st_sample(region_sf, round(sample_size*ss_share[3],0), type="random")
  bpp3<-st_coordinates(bpp_sf3)

  bpp<-rbind(bpp1, bpp2, bpp3)
  bpp$class<-NA
  bpp$class[1:nrow(bpp1)]<-"clustered"
  bpp$class[(nrow(bpp1)+1):(nrow(bpp1)+nrow(bpp2))]<-"regular"
  bpp$class[(nrow(bpp1)+nrow(bpp2)+1):nrow(bpp)]<-"random"

  bpp_sf<-st_as_sf(bpp, coords=c("X", "Y"), crs=st_crs(region_sf), agr="constant")

  plot(bpp_sf, key.pos=1, main="Benford Pattern Spatial Points")
  return(bpp) # poprawić, żeby może na ekranie było jakieś summary/table a nie lista
}
