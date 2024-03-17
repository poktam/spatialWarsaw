#########################
### SpatBenfordTest() ###
#########################
#
# SPRAWDZIĆ CZY NIE OPRZEĆ NA PAKIECIE benford (nowszy) zamiast benford.analysis
#
#' @title Benford's law spatial point pattern conformity test
#'
#' @description
#' Spatial point patterns can conform to Benford's law if the mutual distances between points conform to Benford's law.
#' This occurs in the case of a mixture of regular, random and clustered point patterns composed in stable proportions.
#' The function calculates the mutual Euclidean distances between all points in the dataset and tests whether
#' this distribution conforms to Benford's law.
#'
#' @details
#' The function returns the Benford's law conformity test of geolocated spatial points using the `benford()` function
#' from the `benford.analysis` package. Euclidean distances between points are calculated using the `dist()` function.
#'
#' @name SpatBenfordTest
#' @param data_sf Object in sf of the data.frame class - in the case of a data.frame object, the first and second columns must contain X and Y coordinates.
#' @param sample_size Number of points to be used in the analysis. It must be less than or equal to the number of points in the dataset.
#' Huge datasets can cause computational problems when testing for Benford's law conformity due to the nxn size matrix of mutual distances.
#' @param var_name Column name with additional numeric variable. If empty, a 2d distribution of Euclidean distances is tested, if given, a 3d distribution is tested.
#'
#' @return `SpatBenfordTest()` returns the test for Benford's law conformity.
#'
#' @examples #To be done!!!
#'
#' @export
SpatBenfordTest<-function(data_sf, sample_size, var_name=NULL){
  #Ew. do sprawdzenia czy warunek st_geometry_type(data_sf,FALSE)=="POINT") nie jest zbyt restrykcyjny
  if((inherits(data_sf,"sf") && st_geometry_type(data_sf,FALSE)=="POINT")) crds<-sf::st_coordinates(data_sf)
  else if(inherits(data_sf,"data.frame",TRUE)==1){
    crds<-data_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
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

  # zbadać var_name - jeśli nie ma lub nie istnieje w danych, to pominąć, jeśli istnieje to wyciągnąć kolumnę
  if (length(var_name)>1) {
    var_name<-var_name[1]
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

#' @rdname SpatBenfordTest
#' @export
spatbenfordtest <- SpatBenfordTest


###############################
### SpatBenfordPattern()    ###
###############################
#
# SPRAWDZIĆ CZY NIE OPRZEĆ NA PAKIECIE benford (nowszy) zamiast benford.analysis
# DODAĆ TESTOWANIE WYNIKU ZA POMOCĄ SPATBENFORDPATTERN()
#
#' @title Spatial point pattern generation following Benford's law
#'
#' @description
#' Spatial point patterns can conform to Benford's law if the mutual distances between points conform to Benford's law.
#' This happens in the case of a mixture of regular, random and clustered point patterns composed in stable proportions.
#' The function generates such a point pattern with labels on the original point pattern of each point.
#'
#' @details
#' A natural point pattern (Kopczewska & Kopczewski, 2022) conforming to Benford's law requires a stable proportion
#' of source point patterns: 81% of clustered shifted pattern, 12% of regular pattern and 7% of random point pattern.
#' These proportions have been incorporated into the `SpatBenfordPattern()` function. A clustered shifted point pattern
#' is generated as X and Y coordinates from a normal distribution centred at 0.25 of the spatial range
#' of longitude and latitude. A clustered point pattern centred within the region makes this mixture non-conforming to Benford's law.
#' This natural point pattern can be tested for Benford's law conformity using `SpatBenfordTest()`.
#'
#' @name SpatBenfordPattern
#' @param region_sf Polygon in the `sf` class that defines the boundary for the points to be generated.
#' @param sample_size Number of points to generate inside the boundary. The default value is 5000 and can be changed.
#' Larger datasets may cause computational problems when testing for Benford's law conformity.
#' `St_sample()`, which is used to generate point patterns, may generate fewer points than specified.
#'
#' @return `SpatBenfordPattern()` returns a `data.frame` class object that contains X and Y coordinates in the first
#' two columns and labels on the original distribution of each point in the third column. It also shows the visualisation
#' of a generated mixture of point patterns.
#'
#' @references
#' Kopczewska K, Kopczewski T (2022) Natural spatial pattern—When mutual socio-geo distances between cities follow Benford’s law.
#' PLoS ONE 17(10): e0276450. https://doi.org/10.1371/journal.pone.0276450
#'
#' @seealso [SpatBenfordPattern()]
#'
#' @examples #To be done!!!
#'
#' @export
SpatBenfordPattern<-function(region_sf, sample_size=5000){
  #sprawdzić czy warunek na typ jest ok?
  if(!inherits(region_sf,"sf")) {
    stop("The class of region_sf must only be 'sf'.")
  } else if(!(st_geometry_type(region_sf,FALSE)=="MULTIPOLYGON" || st_geometry_type(region_sf,FALSE)=="POLYGON")){
    stop("The type of region_sf must only be 'MULTIPOLYGON' or 'POLYGON'.")

  }

  ss_share<-c(0.810, 0.120, 0.070) # % of obs. clustered / regular /random
  bb<-st_bbox(region_sf)

  # clustered skewed distribution
  x<-rnorm(sample_size*ss_share[1], bb[1]+0.25*(bb[3]-bb[1]), 0.08*(bb[3]-bb[1]))
  y<-rnorm(sample_size*ss_share[1], bb[2]+0.75*(bb[4]-bb[2]), 0.08*(bb[4]-bb[2]))
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

#' @rdname SpatBenfordPattern
#' @export
spatbenfordpattern <- SpatBenfordPattern
