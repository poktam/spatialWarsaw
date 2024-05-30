#############
### FLE() ###
#############
#
#
#' @title Focal Local Entropy for rasterised data to analyse the density of points
#'
#' @description
#' The function analyses the density of geolocated points by counting the number of points in raster cells
#' and calculating local entropy using neighbouring cells. The counts in the grid cells are normalised
#' and the entropy is calculated for the discretised variable using fixed intervals. The focal mechanism performs
#' the calculations for each grid cell and its neighbours. The nearest neighbours can be given as a number
#' (an odd value e.g. 3,5,7,9) or a radius.
#'
#' @details
#' The function grids the area and counts the number of observations in each grid cell. The counts of the observations are normalised.
#' The focal mechanism selects the neighbourhood for each cell and creates a vector of normalised counts for the analysed cell
#' and its neighbours. The values are discretised into fixed intervals `(breaks=c(-100, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 100))`
#' and the proportion (percentage of observations) in each interval is calculated. These proportions are used as input to the entropy function.
#' The local entropy is calculated for each grid cell and displayed in the raster.
#'
#' The choice of w or r makes a small difference to the result. For w=9, a cell and all its direct neighbours (first row, contiguity)
#' are considered. When setting r, remember that the distance to the centroids of neighbouring cells is closer
#' for east/west/north/south cells and longer for east-north, east-south, west-north and west-south cells.
#' Too short r can exclude the second group from the neighbours list.
#'
#' @name FLE
#' @param points_sf Object in `sf` or the `data.frame` class - in the case of a `data.frame` object, the first and second columns must contain X and Y coordinates.
#' For `data.frame`, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` object.
#' @param region_sf The boundary of the area where the `points_sf` are located, `sf` class object.
#' @param nrows.raster Number of rows in the raster, default 50.
#' @param ncols.raster Number of columns in the raster, default 50.
#' @param w The size of the square moving window (choose an odd value such as 3,5,7,9 etc.); the w and r parameters are mutually exclusive - specify either w or r.
#' @param r The radius of a circular ring in a scale of geo-coordinates; the w and r parameters are mutually exclusive - specify either w or r.
#'
#' @return `FLE()` returns the `terra` class object and its visualisation.
#'
#' @examples
#' my.fle<-FLE(firms_sf, region_sf, nrows.raster=50, ncols.raster=50, w=9)
#' my.fle
#'
#' @export
FLE<-function(points_sf, region_sf, nrows.raster=50, ncols.raster=50, w, r){
  #sprawdzić czy warunek na typ jest ok?
  if(!inherits(region_sf,"sf")) {
    stop("The class of region_sf must only be 'sf'.")
  } else if(!(st_geometry_type(region_sf,FALSE)=="MULTIPOLYGON" || st_geometry_type(region_sf,FALSE)=="POLYGON")){
    stop("The type of region_sf must only be 'MULTIPOLYGON' or 'POLYGON'.")
  }

  #sprawdzić czy wszystko jest potrzebne do tej funkcji
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
      cat("The coordinates for point_sf object of class type data.frame have been tranformed to a geographic coordinate
      system / projection that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".\n",sep="")
    }
  } else if(inherits(points_sf,"data.frame",TRUE)==1){
    crds<-points_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    crds_sf<-st_as_sf(crds,coords = c("X_coord","Y_coord"), crs=st_crs(region_sf), agr="constant")
    cat("The coordinates from the point_sf object of class type data.frame have been assigned a geographic coordinate
    system / projection that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".\n",sep="")
  } else {
    stop("The class of points_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # tylko jeden z dwóch można podać, jeden wymagany
  switch(
    check_exclusive(w, r),
    w = message("`w` was supplied."),
    r = message("`r` was supplied.")
  )

  # czy dawać jeszcze jakieś dodatkowe warunki na nrows.raster, ncols.raster i może jakieś na w lub r? (DO DECYZJI)

  bb<-st_bbox(region_sf) 		# bounding box
  rst<-rast(nrows=nrows.raster, ncols=ncols.raster, xmin=bb[1], ymin=bb[2], xmax=bb[3], ymax=bb[4])	# z terra::
  crds$ones<-rep(1, times=nrow(crds))	# wektor jedynek

  # rastrowanie zmiennej
  rst.var<-rasterize(as.matrix(crds[,1:2]), rst, value=crds$ones,  fun=sum) # z terra:: # czy tu nie powinno być values=?
  rst.var[is.na(rst.var)]<-0
  rst.var<-scale(rst.var)

  # funkcja licząca entropię
  # to są przedziały do zmiennej standaryzowanej
  breaks=c(-100, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 100)

  ent_fun<-function(y) {pi=table(cut(y, breaks=breaks))/sum(table(cut(y, breaks=breaks))); -sum(pi[pi>0]*log(pi[pi>0]))}

  # I wersja z w!!!
  # focal local entropy – w raster
  # wersja korzystająca z parametru w – squared window
  if (!missing(w)) {
    focal_result<-focal(rst.var, w, fun=ent_fun)
    plot(focal_result, main=paste0("Focal local entropy, w=", w))
    plot(st_geometry(region_sf), add=TRUE)
    return (focal_result)
  } else {
    # II wersja z r!!!
    # wersja korzystająca z parametru r – radial window
    w_value <- focalMat(rst.var, r, "circle") # z terra::
    w_value[w_value > 0] <- 1 # replacing weights by 1 to get total
    focal_result<- focal(rst.var, w=w_value, fun = ent_fun)  # z terra::

    plot(focal_result, main=paste0("Focal local entropy, radial window r=", r))
    plot(st_geometry(region_sf), add=TRUE)
    return(focal_result)
  }
}

#' @rdname FLE
#' @export
fle <- FLE
