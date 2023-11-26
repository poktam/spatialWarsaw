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
  my.benford<-benford(as.vector(my.dist))
  plot(my.benford)
  return (my.benford)
}




###############################
### SpatBenfordPattern() ###
###############################
# R code for SpatBenfordPattern()
# generating spatial Benford-like 2D distribution
################################
#
# DODAĆ!!!!!
#
