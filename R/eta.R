#############
### ETA() ###
#############
#
#
#' @title Entropy and Tesselation for Agglomeration (ETA)
#'
#' @description
#' For geolocated point patterns, the degree of agglomeration is calculated. It is expressed as the relative Shannon
#' entropy of the relative areas of tesselation tiles derived as Voronoi polygons. ETA (Entropy and Tesselation
#' for Agglomeration) takes values between 0 and 1. Values close to 1 reflect highly dispersed geolocated points
#' and spatially uniform distribution. The lower the value of ETA, the stronger the agglomeration.
#'
#' @details
#' The function calculates the relative areas of Voronoi polygons (ri) that sum to 1 such that sum(ri)=1.
#' The relative areas (ri) are an input to the Shannon entropy function H, where H=-sum(ri*log(ri)).
#' The maximum entropy Hmax is derived as Hmax=-1/n for n tiles. The relative entropy is calculated as H/Hmax
#' and interpreted as ETA.
#'
#' In the case of spatially uniform distribution (regular location of points), ETA=1. The more points are
#' agglomerated, the lower the ETA. Points agglomerated in a single location should have ETA~0.
#'
#' A number of points may be reduced at the tesselation stage if duplicate locations are identified in the dataset.
#'
#' @name ETA
#' @param points_sf Object in `sf` or the `data.frame` class - in the case of a `data.frame` object, the first and second columns must contain X and Y coordinates.
#' For `data.frame`, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` object.
#' @param region_sf Polygon in the sf class that defines the boundary for points_sf.
#' @param sample_size Sample size, must be less than or equal to the number of points in the dataset (`points_sf` parameter).
#' If `sample_size` is equal to points_sf size, all points are used to construct Voronoi polygons.
#' If `sample_size` is less than the sample size, the functions will perform random sampling.
#' If `sample_size` is greater, it is automatically set to the number of points in the dataset.
#' The number of points affects the speed of Voronoi polygon computation.
#'
#' @return `ETA()` returns a list object and a plot. The list contains three elements:
#' \item{S_ent}{Empirical entropy of the tessellated point pattern.}
#' \item{H_ent}{ETA, degree of agglomeration.}
#' \item{n_points}{Number of points included in calculations.}
#' The plot illustrates the points used in the analysis and the boundary region,
#' Voronoi tesselation and displays values from the list.
#'
#' @references
#' Kopczewska K., 2021, Entropy as measure of agglomeration. Interactions of business locations and housing transactions
#' in Warsaw metropolitan area. \[In] Handbook on "Entropy, Complexity, and Spatial Dynamics: The Rebirth of Theory?",
#' Edward Elgar (editors: Aura Reggiani, Laurie Schintler, Danny Czamanski, Roberto Patuelli)
#'
#' @examples
#' # `sf` package required
#' library(sf)
#'
#' # for random point pattern
#' eta.random<-st_sample(region_sf, 500, type="random")
#' eta.random<-st_as_sf(eta.random)
#' ETA(eta.random, region_sf, 500)
#'
#' # for regular point pattern
#' eta.regular<-st_sample(region_sf, 500, type="regular")
#' eta.regular<-st_as_sf(eta.regular, crs=4326)
#' ETA(eta.regular, region_sf, 500)
#'
#' # for clustered point pattern
#' eta.clust<-data.frame(X=rnorm(500, mean=22.8, sd=0.25), Y=rnorm(500, mean=51.2, sd=0.25))
#' plot(st_geometry(region_sf))
#' points(eta.clust)
#' ETA(eta.clust, region_sf, 500)
#'
#' @export
ETA<-function(points_sf, region_sf, sample_size){
  #sprawdzić czy warunek na typ jest ok?
  if(!inherits(region_sf,"sf")) {
    stop("The class of region_sf must only be 'sf'.\n")
  } else if(!(st_geometry_type(region_sf,FALSE)=="MULTIPOLYGON" || st_geometry_type(region_sf,FALSE)=="POLYGON")){
    stop("The type of region_sf must only be 'MULTIPOLYGON' or 'POLYGON'.")
  }

  # w przypadku gdy oba obiekty są typu sf uzgodnić ich system współrzędnych / projekcję(!!!)
  #Ew. do sprawdzenia czy warunek st_geometry_type(points_sf,FALSE)=="POINT") nie jest zbyt restrykcyjny
  if((inherits(points_sf,"sf") && st_geometry_type(points_sf,FALSE)=="POINT")) {
    # to można uprościć (trochę niepotrzebne wyjęcie współrzędnych i ich przerobienie ponownie), ale na razie zostawimy
    crds<-as.data.frame(st_coordinates(points_sf))
    colnames(crds)<-c("X_coord","Y_coord")
    crds_sf<-st_as_sf(crds,coords = c("X_coord","Y_coord"), crs=st_crs(points_sf), agr="constant")
    if (st_crs(points_sf)!=st_crs(region_sf)) {
      # sprawdzić czy działa przekształcenie
      crds_sf<-st_transform(crds_sf,crs=st_crs(region_sf))
      cat("The coordinates for point_sf object of class type data.frame have been tranformed to a geographic coordinate",
          "system / projection that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".\n",sep="")
    }
  } else if(inherits(points_sf,"data.frame",TRUE)==1){
    crds<-points_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    crds_sf<-st_as_sf(crds,coords = c("X_coord","Y_coord"), crs=st_crs(region_sf), agr="constant")
    cat("The coordinates from the point_sf object of class type data.frame have been assigned ",
      "a geographic coordinate system / projection that matches the projection of the region_sf object: EPSG:", st_crs(region_sf)$epsg,".\n",sep="")
  } else {
    stop("The class of points_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
    }

  # zbadać sample_size i ustawić ew. na wielkość zbioru danych (może potem zmienić, żeby )
  sample_size <- as.integer(sample_size)
  if (sample_size > nrow(points_sf) || sample_size<1) {
    sample_size<-nrow(points_sf)
    cat("Wrong sample size. Sample_size set to:",sample_size,"\n",sep="")
  }

  # przekształcenie do EPSG:3857 zrobić może wcześniej (bo tesselacja lepiej działa w tym układzie - spr?)
  # zmniejszyć liczbę przekształceń
  crds_sf<-st_transform(crds_sf,crs=3857)
  region_sf<-st_transform(region_sf,crs=3857)

  # sample do testowania
  crds_sf_s<-crds_sf[sample(nrow(crds), sample_size, replace=FALSE), ]

  # tesselation - poprawić obiekty (może uprościć)
  crds_sfc_s<-st_geometry(crds_sf_s)
  region_sfc<-st_geometry(region_sf)
  crds_sfc_s_union<-st_union(crds_sfc_s)
  tess_result<-st_voronoi(crds_sfc_s_union, region_sfc)
  tess_result<-st_intersection(st_cast(tess_result), st_union(region_sfc))

  tess_area<-st_area(tess_result)
  tess_area_rel<-tess_area/sum(tess_area)
  S_ent<-sum(-1*tess_area_rel*log(tess_area_rel)) # Shannon entropy
  ent.ref<-log(1/length(tess_area))*(-1)
  H_ent<-S_ent/ ent.ref # Relative H entropy

  # plot with points in blue
  par(mar=c(4,4,4,4))
  plot(st_geometry(tess_result), main="Degree of agglomeration \n entropy of tesselated point pattern",
       sub="Relative H entropy =1 for spatially uniform distribution (no agglomeration)")
  plot(crds_sfc_s, add=TRUE, bg="darkblue", pch=21)
  legend("bottomleft", horiz=FALSE,
         c(paste("Shannon entropy=", round(S_ent,2)), paste("Relative H entropy=", round(H_ent,2)), paste("Number of points =", round(length(tess_area),2))),
         cex=0.75, bty="n")

  list(
    S_ent = as.numeric(S_ent),
    H_ent = as.numeric(H_ent),
    n_points = round(length(tess_area),2)
  )

}

#' @rdname ETA
#' @export
eta <- ETA


