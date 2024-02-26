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

  # cluster your data- xy geocoordinates - to get irregular shapes for sampling
  # czy II parametr zostawić sample_size/100? - co przy b. dużych sample_size?
  points_sf$cluster<-(kmeans(crds, sample_size/100))$cluster

  # to avoid coords in the same place
  crds[,1]<-crds[,1]+rnorm(nrow(crds), 0, sd(crds[,1])/1000)
  crds[,2]<-crds[,2]+rnorm(nrow(crds), 0, sd(crds[,2])/1000)

  # selector selects the rows of observations for every iteration
  # and saves in matrix, later it will be used to recover,
  # on which data the best model was estimated
  # strata() from sampling:: samples obs.from groups (clusters)
  # method="srswor" stands for simple random sampling without replacement (in given column)

  selector<-matrix(NA, nrow=sample_size, ncol=iter)
  for(i in 1:iter){
    #  vec<-sample(nrow(points_sf), sample_size, replace=FALSE)
    selector[,i]<-(strata(points_sf, "cluster", size=rep(100, times=sample_size/100), method="srswor"))$ID_unit
  }

  # objects to store results of iterations
  # (length(var_names)-1) można zastąpić length(attr(terms(eq),"term.labels"))

  n<-ifelse(model_type=="SDM", (length(var_names)-1)*2+1, (length(var_names)-1)+1)

  coef_boot<-matrix(NA, nrow=iter, ncol=n)	# for coefficients
  error_boot<-matrix(NA, nrow=iter, ncol=n) #for std errors
  fitted_boot<-matrix(NA, nrow=sample_size, ncol=iter)	# for fitted values
  y_boot<-matrix(NA, nrow=sample_size, ncol=iter) 	# original values of y

  # AIC.ols, AIC.spatial, time.W, time.model, spatial.coef
  quality_boot<-matrix(NA, nrow=iter, ncol=5)
  colnames(quality_boot)<-c("AIC.ols", "AIC.spatial", "time.W", "time.model", "spatial.coef")

  # main loop – estimation of all models (SEM, SAR, SDM) on the same subsets with time measurement

  cat(iter, " ", model_type, " models will be computed. ",
      "It can take a long while.","\n",sep="")

  for(i in 1:iter){
    dane_temp<-points_sf[selector[,i],]	# subset of data for given iteration

    # W matrix
    crds_temp<-as.matrix(crds[selector[,i], ])
    start_time <- Sys.time()
    knnW_temp<-nb2listw(make.sym.nb(knn2nb(knearneigh(as.matrix(crds_temp), k=knn))))
    end_time <- Sys.time()
    time_W<-difftime(end_time, start_time, units="secs")

    # model estimation
    if (model_type=="SEM") {
      start_time <- Sys.time()
      model_temp<-errorsarlm(eq, data=dane_temp, knnW_temp, method="LU")
      end_time <- Sys.time()
      time_model<- difftime(end_time, start_time, units="secs")

    } else if (model_type=="SAR") {
      start_time <- Sys.time()
      model_temp<-lagsarlm(eq, data=dane_temp, knnW_temp, method="LU")
      end_time <- Sys.time()
      time_model<- difftime(end_time, start_time, units="secs")

    } else if (model_type=="SDM") {
      start_time <- Sys.time()
      model_temp<-lagsarlm(eq, data=dane_temp, knnW_temp, method="LU", type="mixed")
      end_time <- Sys.time()
      time_model<- difftime(end_time, start_time, units="secs")
    }

    # saving the results into appropriate objects
    # sprawdzić czy nie wystarczy coef_boot[i,]<-model_temp$coefficients itd.
    coef_boot[i,1:length(model_temp$coefficients)]<-model_temp$coefficients
    error_boot[i, 1:length(model_temp$coefficients)]<-model_temp$rest.se
    fitted_boot[,i]<-model_temp$fitted.values
    y_boot[,i]<-as.matrix(st_drop_geometry(dane_temp[,var_names[1]]))
    quality_boot[i,1]<-model_temp$AIC_lm.model # AIC.ols
    quality_boot[i,2]<-AIC(model_temp) # AIC.spatial
    quality_boot[i,3]<-time_W
    quality_boot[i,4]<-time_model
    quality_boot[i,5]<-ifelse(is.null(model_temp$rho)==TRUE, model_temp$lambda, model_temp$rho)
  }


  # selection of the best model with PAM
  clust_mod<-pam(cbind(coef_boot, quality_boot[,5]),1)

  dane_best<-points_sf[selector[,clust_mod$id.med],] # data which were used in estimation of the best model
  crds_best<-crds[selector[,clust_mod$id.med],]
  RAMSE_best<-(sum((y_boot[ , clust_mod$id.med] - fitted_boot[ , clust_mod$id.med])^2)/sample_size)^(0.5)
  knnW_best<-nb2listw(make.sym.nb(knn2nb(knearneigh(as.matrix(crds_best), k=knn))))

  # best model estimation
  if (model_type=="SEM") {
    model_best<-errorsarlm(eq, data=dane_best, knnW_best, method="LU")

  } else if (model_type=="SAR") {
    model_best<-lagsarlm(eq, data=dane_best, knnW_best, method="LU")

  } else if (model_type=="SDM") {
    model_best<-lagsarlm(eq, data=dane_best, knnW_best, method="LU", type="mixed")
  }

  # SPRAWDZIĆ, czy kolumny mają się tak nazywać
  colnames(coef_boot)<-names(model_best$coefficients)
  colnames(error_boot)<-names(model_best$coefficients)

  list(coef.boot = coef_boot,
       error.boot = error_boot,
       quality.boot = quality_boot,
       dane.best = dane_best,
       model.best = model_best,
       RAMSE.best = RAMSE_best)
  #może przygotować klasę do tego outputu?

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

