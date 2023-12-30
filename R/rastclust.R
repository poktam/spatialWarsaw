######################
### rastClustGWR() ###
######################
#
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name rastClustGWR
#' @param points_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords).
#' When using a simple data.frame, make sure that the coordinates of the points are in the same
#' coordinate system / projection as the `region_sf` object. NOTE! Data must be from a single period!
#' @param eq an object of class [stats::formula()] (or one that can be coerced to that class):
#' a symbolic description of the model to be used.
#' @param region_sf do opisu (obiekt sf ale jako region)
#' @param nrows.raster raster row dimensions, default 50
#' @param ncols.raster raster column dimensions, default 50
#' @param nc number of clusters
#' @param bw do opisu
#' @param adaptive is adaptive (do opisu), default=FALSE
#'
#' @examples #To be done!!!
#'
#' @return `rastClustGWR()` returns ... to be done.
#'
#' @export
rastClustGWR<-function(points_sf, eq, region_sf, nrows.raster=50, ncols.raster=50, nc, bw, adaptive=FALSE){
  #sprawdzić czy warunek na typ jest ok?
  if(!inherits(region_sf,"sf")) {
    stop("The class of region_sf must only be 'sf'.\n")
  } else if(!(st_geometry_type(region_sf,FALSE)=="MULTIPOLYGON" || st_geometry_type(region_sf,FALSE)=="POLYGON")){
    stop("The type of region_sf must only be 'MULTIPOLYGON' or 'POLYGON'.")
  }

  #Ew. do sprawdzenia czy warunek st_geometry_type(points_sf,FALSE)=="POINT") nie jest zbyt restrykcyjny
  if((inherits(points_sf,"sf") && st_geometry_type(points_sf,FALSE)=="POINT")) {
    cat("Points_sf was detected as an object of class sf.\n")
    if (st_crs(points_sf)!=st_crs(region_sf)) {
      # sprawdzić czy działa przekształcenie
      points_sf<-st_transform(points_sf,crs=st_crs(region_sf))
      cat("The coordinates for point_sf object of class type data.frame have been tranformed to a geographic coordinate
      system / projection that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".\n",sep="")
    }
  } else if(inherits(points_sf,"data.frame",TRUE)==1){
    colnames(points_sf)[1:2]<-c("X_coord","Y_coord")
    points_sf<-st_as_sf(points_sf,coords = c("X_coord","Y_coord"), crs=st_crs(region_sf), agr="constant")
    cat("The coordinates from the point_sf object of class type data.frame have been assigned a geographic coordinate system / projection
    that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".
    The points_sf object has been converted to sf class.\n",sep="")
  } else {
    stop("The class of points_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  var_names<-all.vars(eq)
  m <- match(gsub(" ", ".", var_names), colnames(points_sf))
  if (any(is.na(m))) {
    stop("Variable names are incorrect.")
  }

  # czy dawać jeszcze jakieś dodatkowe warunki na nrows.raster, ncols.raster, nc, bw? (DO DECYZJI)

  #GW wymaga sp :(
  crds<-st_coordinates(points_sf)
  points_sp<-as_Spatial(points_sf)



  bb<-st_bbox(region_sf) 		# bounding box
  rst<-rast(nrows=nrows.raster, ncols=ncols.raster, xmin=bb[1], ymin=bb[2], xmax=bb[3], ymax=bb[4])	# z terra::

  dMat<-gw.dist(crds) # czasem długo trwa!

  # jeśli użytkownik podał wartość, to nie liczyć tego – długo trwa
  if (missing(bw) || is.null(bw)) {
    cat("Parameter bw not specified. \n")
    bw<-bw.gwr(eq, data=points_sp, kernel='gaussian', adaptive=adaptive, dMat=dMat)
  }

  # w gwr.basic jest longlat=F- czy nie trzeba dla współrzędnych
  modelGWR<-gwr.basic(eq, data=points_sp, kernel='gaussian', adaptive=adaptive, bw=bw, dMat=dMat)

  # clusters of coefficients – we use k-means
  end<-match("y",names(modelGWR$SDF))-1 # można się pozbyć zmiennej end (ew.)
  km.clust<-kmeans(as.data.frame(modelGWR$SDF[,1:end]), nc)
  rst.var<-rasterize(crds, rst, value=km.clust$clust,  fun=median) #terra

  plot(rst.var, main="Rasters with median ID \n k-means clusters of GWR coefficients", cex.main=0.9, mar=c(2,2,2.5,2))
  plot(st_geometry(region_sf), add=TRUE)
  return (rst.var)
}
