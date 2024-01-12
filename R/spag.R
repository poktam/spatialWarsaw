##############
### SPAG() ###
##############
#
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name SPAG
#' @aliases SPAG spag
#' @param points_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords). When using a simple data.frame, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` object. (OD RAZU DANE Z WYBRANEJ ZMIENNEJ NP. SEKTORA)
#' @param size_var The name of the variable describing the size of the analysed objects (e.g. companies). (SPRAWDZIĆ opis)
#' @param region_sf do opisu (obiekt sf ale jako region)
#' @examples #To be done!!!
#'
#' @return `SPAG()` returns ... to be done.
#'
#' @export
SPAG<-function(points_sf, size_var, region_sf){
  if(!inherits(region_sf,"sf")) {
    stop("The class of region_sf must only be 'sf'.\n")
  } else if(!(st_geometry_type(region_sf,FALSE)=="MULTIPOLYGON" || st_geometry_type(region_sf,FALSE)=="POLYGON")){
    stop("The type of region_sf must only be 'MULTIPOLYGON' or 'POLYGON'.")
  }

  #Ew. do sprawdzenia czy warunek st_geometry_type(points_sf,FALSE)=="POINT") nie jest zbyt restrykcyjny
  if((inherits(points_sf,"sf") && st_geometry_type(points_sf,FALSE)=="POINT")) {
    cat("Points_sf was detected as an object of class sf.\n", sep = "")
    if (st_crs(points_sf)!=st_crs(region_sf)) {
      # sprawdzić czy działa przekształcenie
      points_sf<-st_transform(points_sf,crs=st_crs(region_sf))
      cat("The coordinates for point_sf object of class type data.frame have been tranformed to a geographic coordinate
      system / projection that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".\n",sep="")
    }
  } else if(inherits(points_sf,"data.frame",TRUE)==1){
    colnames(points_sf)[1:2]<-c("X_coord","Y_coord")
    points_sf<-st_as_sf(points_sf,coords = c("X_coord","Y_coord"), crs=st_crs(region_sf), agr="constant")
    cat("The coordinates from the point_sf object of class type data.frame have been assigned ",
        "a geographic coordinate system / projection that matches the projection of the region_sf object: EPSG:", st_crs(region_sf)$epsg,".\n",
        "The points_sf object has been converted to sf class.\n",sep="")
  } else {
    stop("The class of points_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # SPRAWDZIĆ CZY NIE LEPIEJ OD RAZU PRZEKSZTAŁCIĆ REGION I PUNKTY_SF DO EPSG:4326 - CHYBA W TERRA:: LEPIEJ TO DZIAŁA


  # zbadać size_var - jeśli nie ma lub nie istnieje w danych, to pominąć, jeśli istnieje to wyciągnąć kolumnę
  if (length(size_var)>1) {
    size_var<-size_var[1]
    cat("Parameter size_var longer than 1. The first element has been selected: ",size_var,"\n",sep="")
  }

  m <- match(gsub(" ", ".", size_var), colnames(points_sf))

  if (is.na(m) || is.null(size_var)) {
    stop("Unknown variable name or variable name not specified.")
  }

  points_sf_s<-points_sf[,m]
  colnames(points_sf_s)[1]<-"size"
  points_sf_s$size<-round(points_sf_s$size)


  # lekki porządek poniżej do zrobienia (choć działa)
  sizes<-unique(round(as.data.frame(points_sf[,m])[,1]))
  firm_num<-as.data.frame(table(round(as.data.frame(points_sf[,m])[,1])))
  colnames(firm_num)[1]<-"size"
  firm_num$rel_size<-sizes/min(sizes)
  firm_num$r_multipl<-sqrt(firm_num$rel_size)

  area_region<-st_area(region_sf) #denominator of nomin.coverage.sel

  f.char<-character()
  for(i in 1:nrow(firm_num)){
    f.char<-paste(f.char,paste("firm_num[",i,",2]*pi*firm_num[",i,",3]*r^2",sep=''),sep='+')
  }
  eval(parse(text = paste('f.all <- function(r) { return(' , f.char , '-as.numeric(area_region))}', sep='')))

  res_root<-uniroot(f.all, c(0, 100000), tol = 1e-35)

  points_sf_s<-merge(points_sf_s,firm_num[,c(1,4)],by="size")
  points_sf_s$rad<-points_sf_s$r_multipl*res_root$root
  points_sf_s$area<-pi*(points_sf_s$rad)^2

  points_vec_s<-vect(points_sf_s)
  circles_vec_s<-buffer(points_vec_s, width=points_vec_s$rad) # single circles
  circles_vec_s<-aggregate(circles_vec_s) #merged circles

  # random matrix of distances
  # MOŻE 3000 UZALEŻNIĆ OD WIELKOŚCI OBSZARU - JAKOŚ PROPORCJONALNIE
  selector_rnd<-sample(1:nrow(points_sf), min(nrow(points_sf_s),3000), replace = FALSE)
  points_sf_rnd_crds<-st_coordinates(points_sf_s[selector_s,]) # any random points from the sample
  distance_rnd<-dist(points_sf_rnd_crds)

  # counters for analysed sector
  counter_distance<-mean(distance_rnd)
  counter_coverage<-sum(points_sf_s$area)
  counter_overlap<-expanse(circles_vec_s)

  # theoretical part
  # equal location of selected points
  # we use object ‘punkty’ with uniform distribution
  # number of points is given wi ‘nn.real’


  points_sf_theo<-st_sample(region_sf, 5000, type="regular")
  st_crs(points_sf_theo)<-st_crs(region_sf)
  # points_sf_theo<-st_transform(points_sf_theo, 4326) 	# NA RAZIE BEZ KONWERSJI - spr. konwersja do WGS84 (bo lepiej działa)

  theo_rad<-sqrt(area_region/(length(points_sf_theo)*pi))

  # distances for selected theoretical locations   random selection
  selector_theo<-sample(1:length(points_sf_theo), min(100, length(points_sf_theo)), replace=FALSE)
  points_sf_theo_rnd_crds<-st_coordinates(points_sf_theo[selector_theo])
  distance_theo<-dist(points_sf_theo_rnd_crds)

  nomin_distance<-mean(distance_theo)
  nomin_coverage<-area_region
  nomin_overlap<-length(points_sf_theo)*pi*theo_rad^2 #in fact area of region

  # CZY nomin_coverage = nomin_overlap = area_region ???? SPRAWDZIĆ???

  #SPAG
  i_coverage<-as.numeric(counter_coverage/nomin_coverage)
  i_distance<-as.numeric(counter_distance/nomin_distance)
  i_overlap<-as.numeric(counter_overlap/nomin_overlap)
  SPAG_res<- i.coverage * i.distance * i.overlap

  # MOŻE PROBLEM Z LEGENDĄ - sprawdzić
  par(mar=c(2,2,2,2)+0.1)
  plot(st_geometry(WOJ_LUB))
  plot(circles_vec_s, col="lightblue", add=TRUE)
  legend("bottomleft", legend=c(paste("i.coverage=",round(i_coverage,2)),
                                paste("i.distance=",round(i_distance,2)),
                                paste("i.overlap=",round(i_overlap,2)),
                                paste("SPAG=",round(SPAG_res,3)),
                                paste("n obs.=",round(nrow(points_sf_s),2))), cex=0.85, bty="n")

  title(main="SPAG measure for the analysed data", cex.main=0.9)

  list(i.coverage = i_coverage,
       i.distance = i_distance,
       i.overlap = i_overlap,
       SPAG = SPAG_res,
       n.obs = nrow(points_sf_s))
}
#' @export
spag <- SPAG


