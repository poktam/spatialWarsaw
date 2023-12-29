#############
### FLE() ###
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
#' @name FLE
#' @aliases FLE fle
#' @param points_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords). When using a simple data.frame, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` object.
#' @param region_sf do opisu (obiekt sf ale jako region)
#' @param nrows.raster raster row dimensions, default 50
#' @param ncols.raster raster column dimensions, default 50
#' @param w size of squared moving window (choose 3,5,7,9 etc), w and r parameters are mutually exclusive
#' @param r radius of circular ring, w and r parameters are mutually exclusive
#'
#' @examples #To be done!!!
#'
#' @return `FLE()` returns ... to be done.
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

  # tylko jeden z dwóch można podać (sprawdzić) - może jeszcze, że oba nullami nie może być
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
  rast.var<-rasterize(as.matrix(crds[,1:2]), rst, value=crds$ones,  fun=sum) # z terra::
  rast.var.vec<-as.vector(rast.var$sum)
  rast.var.vec[is.na(rast.var.vec)==TRUE]<-0

  rst2<-rst						# kopia obiektu (po co?)
  terra::values(rst2)<-scale(rast.var.vec)		# normalizacja

  # funkcja licząca entropię
  # to są przedziały do zmiennej standaryzowanej, zawsze działają
  breaks=c(-100, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 100)

  ff<-function(y) {pi=table(cut(y, breaks=breaks))/sum(table(cut(y, breaks=breaks))); -sum(pi[pi>0]*log(pi[pi>0]))}

  # I wersja z w!!!
  # focal local entropy – w raster
  # wersja korzystająca z parametru w – squared window
  if (!missing(w)) {
    gg<-focal(rst2, w, fun=ff)
    plot(gg, main=paste0("Focal local entropy, w=", w))
    plot(st_geometry(region_sf), add=TRUE)
    return (gg)
  } else {
    # II wersja z r!!!
    # wersja korzystająca z parametru r – radial window
    w_terra <- focalMat(rst2, r, "circle") # z terra::
    w_terra[w_terra > 0] <- 1 # replacing weights by 1 to get total
    focal_terra<- focal(rst2, w=w_terra, fun = ff)  # z terra::

    plot(focal_terra, main=paste0("Focal local entropy, radial window r=", r))
    plot(st_geometry(region_sf), add=TRUE)
    return(focal_terra)
  }
}
#' @export
fle <- FLE
