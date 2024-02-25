#####################
### BootSpatReg() ###
#####################
#
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name BootSpatReg
#' @aliases BootSpatReg
#' @param points_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords).
#' @param iter Number of iterations
#' @param sample_size Sample size, must be less than or equal to the number of points in the dataset (`points_sf` parameter).
#' If `sample_size` is larger, it is automatically set to the number of points in the dataset.
#' We suggest that a value greater than 800 is not used for reasons of computational efficiency. (SPRAWDZIĆ)
#' @param eq an object of class [stats::formula()] (or one that can be coerced to that class):
#' a symbolic description of the model to be used.
#' @param model_type one of: "SAR","SDM","SEM". Default set to "SDM"
#' @param knn chosen knn
#'
#' @examples #To be done!!!
#'
#' @return `BootSpatReg()` returns ... to be done.
#'
#' @export
BootSpatReg<-function(points_sf, iter, sample_size, eq, model_type="SDM", knn){
  #Ew. do sprawdzenia czy warunek st_geometry_type(data_sf,FALSE)=="POINT") nie jest zbyt restrykcyjny
  if((inherits(points_sf,"sf") && st_geometry_type(points_sf,FALSE)=="POINT")) {
    crds<-as.data.frame(st_coordinates(points_sf))
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Points_sf was detected as an object of class sf.\n", sep = "")
  }  else if(inherits(points_sf,"data.frame",TRUE)==1){
    crds<-points_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Points_sf was detected as an object of class data.frame.\n", sep = "")
  }  else {
    stop("The class of data_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # zbadać sample_size i ustawić ew. na wielkość zbioru danych (może potem zmienić, żeby )
  sample_size <- as.integer(sample_size)
  if (sample_size > nrow(points_sf) || sample_size<1) {
    sample_size<-nrow(points_sf)
    cat("Wrong sample size. Sample_size set to:",sample_size,"\n",sep="")
  }

  var_names<-all.vars(eq)
  m <- match(gsub(" ", ".", var_names), colnames(points_sf))
  if (any(is.na(m))) {
    stop("Variable names in eq are incorrect.")
  }

  # zbadać model_type - jeśli nie ma lub nie istnieje w danych, to pominąć, jeśli istnieje to wyciągnąć kolumnę
  if (length(model_type)>1) {
    model_type<-model_type[1]
    cat("Parameter model_type longer than 1. The first element has been selected: ",model_type,"\n",sep="")
  }

  if (!(model_type %in% c("SAR","SDM","SEM"))) {
    stop("Unknown model type. Must be one of: SAR, SDM, SEM.")
  }

  # knn - warunek
  if(!(is.numeric(knn))) {
    stop("knn is to a numerical value.")
  } else if (length(knn)>1) {
    knn<-knn[1]
    cat("knn should be a single value. The first one given was chosen.\n", sep = "")
  }
  knn<-round(knn)

  #iter - warunek
  if(iter<0){
    iter<-100
    cat("iter should be >0. Set to default value 100.\n", sep = "")
  }
  iter<-round(iter)










  # MOŻE PROBLEM Z LEGENDĄ - sprawdzić
  par(mar=c(2,2,2,2)+0.1)
  plot(st_geometry(region_sf))
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

#######################
### ApproxSERoot2() ###
#######################
#
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name ApproxSERoot2
#' @aliases ApproxSERoot2
#' @param points_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords). When using a simple data.frame, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` object. (OD RAZU DANE Z WYBRANEJ ZMIENNEJ NP. SEKTORA)
#' @param size_var The name of the variable describing the size of the analysed objects (e.g. companies). (SPRAWDZIĆ opis)
#' @param region_sf do opisu (obiekt sf ale jako region)
#' @examples #To be done!!!
#'
#' @return `ApproxSERoot2()` returns ... to be done.
#'
#' @export
ApproxSERoot2<-function(points_sf, size_var, region_sf){
}

#######################
### SpatPredTess() ###
#######################
#
#
#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description
#'
#' Linijka nr 3 - details
#'
#'
#' @name SpatPredTess
#' @aliases SpatPredTess
#' @param points_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords). When using a simple data.frame, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` object. (OD RAZU DANE Z WYBRANEJ ZMIENNEJ NP. SEKTORA)
#' @param size_var The name of the variable describing the size of the analysed objects (e.g. companies). (SPRAWDZIĆ opis)
#' @param region_sf do opisu (obiekt sf ale jako region)
#' @examples #To be done!!!
#'
#' @return `SpatPredTess()` returns ... to be done.
#'
#' @export
SpatPredTess<-function(points_sf, size_var, region_sf){
}

