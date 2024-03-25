#####################
### BootSpatReg() ###
#####################
#
#' @title Spatial bootstrapped regression model for handling large data sets
#'
#' @description
#' The function uses the bootstrap approach to estimate a set of spatial models, each on a different subset of the data - the function selects the rows
#' and the columns are fixed according to the equation. The best model is selected using the Partitioning Around Medoid (PAM) algorithm for the k=1 cluster -
#' it is the most central model in a multidimensional setting. The function reports which observations are associated with the best model.
#'
#' @details
#' For each iteration specified in the iter argument, the function selects the random subset of observations, estimates the model and stores the results.
#' The best model is selected as the most central using multidimensional medoid for all estimated coefficients. PAM searches for the model that is closest
#' to the other models on all coefficients.
#'
#' The models are estimated more quickly on subsamples than on a full sample. As shown in Kopczewska (2023), coefficients in subsample models are consistent -
#' there is no significant difference between coefficients from full and subsample models. The standard errors of the coefficients in a subsample model
#' are generally higher than in a full sample model. However, they can be approximated by the √2 rule - the standard error in a model based on a sample twice as large
#' is √2 as small. This approximation can be easily checked using the [ApproxSERoot2()] function. This property allows estimating spatial econometric models
#' for large data using random subsamples without losing precision.
#'
#' `BootSpatReg()` is a time-efficient, methodologically correct approach to deal with the non-scalability of the spatial weight matrix W.
#'
#' @name BootSpatReg
#' @param points_sf do opisu (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords).
#' @param iter The number of iterations.
#' @param sample_size The sample size, must be less than or equal to the number of points in the dataset (`points_sf` parameter). If `sample_size` is greater, it is automatically set
#' to the number of points in the dataset.For the first trial, try a value around 1000 for computational efficiency.
#' @param eq An object of class [stats::formula()] (or one that can be coerced to that class) that defines the equation for the model: a symbolic description of the model to be used.
#' @param model_type The type of spatial econometric model, one of "SAR", "SDM", "SEM". *Soon to be updated for use with "SDEM", "SAC" models.*
#' @param knn The number of k nearest neighbours (knn) used to construct the spatial weight matrix based on the k nearest neighbours criterion.
#'
#' @return `BootSpatReg()` returns the following list object:
#' \item{coef.boot}{A `data.frame` with the coefficients of all iterations.}
#' \item{error.boot}{A `data.frame` with the standard errors of the coefficients from all iterations.}
#' \item{quality.boot}{A set of quality metrics for each iteration: Akaike Information Criterion of OLS model (`AIC.ols`) and selected spatial model (`AIC.spatial`),
#' time of computation of spatial weight matrix (`time.W`) and estimation of spatial model (`time.model`), selected spatial coefficient, rho or lambda (`spatial.coef`).}
#' \item{data.best}{A subset of the data used to estimate the best model, reported as an `sf` object; the rownames (`outome$data.best`) command allows you to obtain the IDs of the rows
#' used in the iteration that produced the best model.}
#' \item{knnW.best}{A spatial weight matrix used to estimate the best spatial model, constructed using the k nearest neighbours criterion for k specified by the user as input.}
#' \item{model.best}{An object of the best spatial model, selected with PAM from the set of bootstrap models.}
#' \item{RAMSE.best}{Root Mean Square Error (RAMSE) of the best spatial model.}
#'
#' @references
#' Kopczewska, K. (2023). Spatial bootstrapped microeconometrics: Forecasting for out‐of‐sample geo‐locations in big data.
#' Scandinavian Journal of Statistics.
#'
#' @examples #To be done!!!
#'
#' @export
BootSpatReg<-function(points_sf, iter, sample_size, eq, model_type, knn){
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

  model_best$call$formula<-eq # wpisanie formuły do modelu

  # SPRAWDZIĆ, czy kolumny mają się tak nazywać
  colnames(coef_boot)<-names(model_best$coefficients)
  colnames(error_boot)<-names(model_best$coefficients)

  list(coef.boot = coef_boot,
       error.boot = error_boot,
       quality.boot = quality_boot,
       data.best = dane_best,
       knnW.best = knnW_best,
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
#' @param model_spatial Spatial model object to analyse
#' @examples #To be done!!!
#'
#' @return `ApproxSERoot2()` returns ... to be done.
#'
#' @export
ApproxSERoot2<-function(model_spatial){
  if(!(inherits(model_spatial,"Sarlm"))) {
    stop("The class of model_spatial must be 'Sarlm' only.")
  }

  #można rozszerzyć na wybraną w argumentach funkcji liczbę potęg 2.
  mp<-length(model_spatial$y)*c(1, 2, 4, 8, 16, 32) # Sample size multiplied by successive powers of 2
  result<-matrix(NA, nrow=length(model_spatial$coefficients), ncol=7)
  colnames(result)<-c("coef", paste0("SE_",mp[1]), paste0("SE_",mp[2]), paste0("SE_",mp[3]), paste0("SE_",mp[4]), paste0("SE_",mp[5]), paste0("SE_",mp[6]))
  rownames(result)<-names(model_spatial$coefficients)
  result[,1]<-model_spatial$coefficients

  # wykorzystać!!! (bez pętli)
  # as.matrix(model_spatial$rest.se)%*%t(as.matrix(c(1, 1/(2^0.5)^1, 1/(2^0.5)^2, 1/(2^0.5)^3, 1/(2^0.5)^4, 1/(2^0.5)^5)))
  for(i in 1:length(model_spatial$coefficients)){
    er<-model_spatial$rest.se[i]
    result[i, 2:7]<-c(er, er/(2^0.5)^1, er/(2^0.5)^2, er/(2^0.5)^3, er/(2^0.5)^4, er/(2^0.5)^5)
  }

  # plot (ew. doszlifować)
  plot(mp, result[1,2:7], type="l", xlab="size of dataset", ylab="SE of beta coefficient",
       ylim=c(min(result[,2:7]),max(result[,2:7])), col="red", lty=1, lwd=2)
  text(mp[1]/2,result[1,2], "1", cex=0.8)

  for(i in 2:length(model_spatial$coefficients)){
    lines(mp, result[i, 2:7], lty=1, col="red", lwd=1)
    text(mp[1]/2,result[i,2], i, cex=0.8)}

  title(main="Extrapolated SE of model coefficients based on sqrt(2) rule")
  legend("topright", paste0(1:length(model_spatial$coefficients), ":", names(model_spatial$coefficients) ), bty="n", cex=0.8)

  result
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
#' @param model_spatial Spatial model
#' @param points_spatial_sf Data on which the model_spatial was estimated (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords).
#' When using a simple data.frame, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` object.
#' @param knnW Spatial weighting matrix (`listw` class) used to estimate model_spatial
#' @param points_new_sf New data for prediction (all used for a forecast). NOTE: The new data must have the same class and structure as
#' the `points_spatial_sf` object. (obiekt sf lub data.frame - w data frame 1 kolumna musi być X coords, druga kolumna Y coords).
#' When using a simple data.frame, make sure that the coordinates of the points are in the same coordinate system / projection as the `region_sf` and `points_spatial_sf` objects.
#' @param region_sf do opisu (obiekt sf ale jako region)
#' @examples #To be done!!!
#'
#' @return `SpatPredTess()` returns ... to be done.
#'
#' @export
SpatPredTess<-function(model_spatial, points_spatial_sf, knnW, points_new_sf, region_sf){
  if(!(inherits(model_spatial,"Sarlm"))) {
    stop("The class of model_spatial must be 'Sarlm' only.")
  }

  #sprawdzić czy warunek na typ jest ok?
  if(!inherits(region_sf,"sf")) {
    stop("The class of region_sf must only be 'sf'.\n")
  } else if(!(st_geometry_type(region_sf,FALSE)=="MULTIPOLYGON" || st_geometry_type(region_sf,FALSE)=="POLYGON")){
    stop("The type of region_sf must only be 'MULTIPOLYGON' or 'POLYGON'.")
  }

  # w przypadku gdy oba obiekty są typu sf uzgodnić ich system współrzędnych / projekcję(!!!)
  #Ew. do sprawdzenia czy warunek st_geometry_type(points_spatial_sf,FALSE)=="POINT") nie jest zbyt restrykcyjny
  if((inherits(points_spatial_sf,"sf") && st_geometry_type(points_spatial_sf,FALSE)=="POINT")) {
    # to można uprościć (trochę niepotrzebne wyjęcie współrzędnych i ich przerobienie ponownie), ale na razie zostawimy
    if (st_crs(points_spatial_sf)!=st_crs(region_sf)) {
      # sprawdzić czy działa przekształcenie
      points_spatial_sf<-st_transform(crds_spatial_sf,crs=st_crs(region_sf))
      cat("The coordinates for point_spatial_sf object of class type data.frame have been tranformed to a geographic coordinate",
          "system / projection that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".\n",sep="")
    }
    crds_spatial<-as.data.frame(st_coordinates(points_spatial_sf))
    colnames(crds_spatial)<-c("X_coord","Y_coord")
    crds_spatial_sf<-st_as_sf(crds_spatial,coords = c("X_coord","Y_coord"), crs=st_crs(points_spatial_sf), agr="constant")
  } else if(inherits(points_spatial_sf,"data.frame",TRUE)==1){
    points_spatial_sf<-st_as_sf(points_spatial_sf,coords = c(1,2), crs=st_crs(region_sf), agr="constant")
    crds_spatial<-as.data.frame(st_coordinates(points_spatial_sf))
    colnames(crds_spatial)<-c("X_coord","Y_coord")
    crds_spatial_sf<-st_as_sf(crds_spatial,coords = c("X_coord","Y_coord"), crs=st_crs(points_spatial_sf), agr="constant")
    cat("The coordinates from the point_spatial_sf object of class type data.frame have been assigned ",
        "a geographic coordinate system / projection that matches the projection of the region_sf object: EPSG:", st_crs(region_sf)$epsg,".\n",sep="")
  } else {
    stop("The class of points_spatial_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  if((inherits(points_new_sf,"sf") && st_geometry_type(points_new_sf,FALSE)=="POINT")) {
    # to można uprościć (trochę niepotrzebne wyjęcie współrzędnych i ich przerobienie ponownie), ale na razie zostawimy
    if (st_crs(points_new_sf)!=st_crs(region_sf)) {
      # sprawdzić czy działa przekształcenie
      points_new_sf<-st_transform(crds_spatial_sf,crs=st_crs(region_sf))
      cat("The coordinates for point_new_sf object of class type data.frame have been tranformed to a geographic coordinate",
          "system / projection that matches the projection of the region_sf object: EPSG:",st_crs(region_sf)$epsg,".\n",sep="")
    }
    crds_new<-as.data.frame(st_coordinates(points_new_sf))
    colnames(crds_new)<-c("X_coord","Y_coord")
    crds_new_sf<-st_as_sf(crds_new,coords = c("X_coord","Y_coord"), crs=st_crs(points_new_sf), agr="constant")
  } else if(inherits(points_new_sf,"data.frame",TRUE)==1){
    points_new_sf<-st_as_sf(points_new_sf,coords = c(1,2), crs=st_crs(region_sf), agr="constant")
    crds_new<-as.data.frame(st_coordinates(points_new_sf))
    colnames(crds_new)<-c("X_coord","Y_coord")
    crds_new_sf<-st_as_sf(crds_new,coords = c("X_coord","Y_coord"), crs=st_crs(points_new_sf), agr="constant")
    cat("The coordinates from the point_spatial_sf object of class type data.frame have been assigned ",
        "a geographic coordinate system / projection that matches the projection of the region_sf object: EPSG:", st_crs(region_sf)$epsg,".\n",sep="")
  } else {
    stop("The class of points_spatial_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  if(!(inherits(knnW,"listw"))) {
    stop("The class of knnW must be 'listw'.")
  }

  # SPRAWDZIĆ/DODAĆ sprawdzenie, czy struktura (nazwy kolumn) points_spatial_sf  odpowiadają points spatial new_sf.
  if (!identical(colnames(points_spatial_sf),colnames(points_new_sf))) {
    stop("Names of variables in points_spatial_sf and points_new_sf are not the same.")
  }

  size_new <- nrow(points_new_sf)

  # For tesselation only
  crds_tess_sf<-st_transform(crds_spatial_sf,crs=3857)
  crds_tess_new_sf<-st_transform(crds_new_sf,crs=3857)
  region_tess_sf<-st_transform(region_sf,crs=3857)

  # tesselation (może uprościć)
  crds_tess_sfc<-st_geometry(crds_tess_sf)
  region_tess_sfc<-st_geometry(region_tess_sf)
  crds_tess_sfc_union<-st_union(crds_tess_sfc)
  tess_result<-st_voronoi(crds_tess_sfc_union, region_tess_sfc)
  tess_result<-st_intersection(st_cast(tess_result), st_union(region_tess_sfc))

  tess_result_sf<-st_sf(tess_result)
  ppi<-st_intersects(crds_tess_new_sf, tess_result_sf) # prediction points indicator

  var_dep<-all.vars(model_spatial$call$formula)[1]

  forecasts_result<-matrix(NA, nrow=size_new, ncol=5)
  colnames(forecasts_result)<-c("predicted y", "real y", "crds x", "crds y", "(predY-realY)^2")

  # loop for forecasts for new points - separate match for each point
  for(i in 1:size_new){
    points_pred_sf<-points_spatial_sf
    points_pred_sf[unlist(ppi[i]),] <- points_new_sf[i, ]
    rownames(points_pred_sf)<-1:nrow(points_pred_sf)

    # prediction for out-of-sample calibrated selected model
    pred_temp<-predict(model_spatial, newdata=points_pred_sf, listw=knnW, legacy.mixed=TRUE)
    forecasts_result[i,1]<- pred_temp[unlist(ppi[i])] # predicted y for a new point
    forecasts_result[i,2]<- st_drop_geometry(points_new_sf[i, var_dep])[1,1]		# empirical y for a new point
    forecasts_result[i,3:4]<- st_coordinates(points_new_sf[i, var_dep]) 		# x,y coordinates
  }

  forecasts_result[,5]<-(forecasts_result[,1]-forecasts_result[,2])^2
  RAMSE_pred<-(mean(forecasts_result[,5]))^0.5

  # plotting tessellation + points for predictions

  par(mar=c(1,1,2,1))
  plot(st_transform(tess_result,crs=st_crs(region_sf)), main="Points used in model and for prediction")	# tessellation plot
  plot(st_geometry(crds_spatial_sf), bg="blue", pch=21, cex=0.5, add=TRUE)	# points from model
  plot(st_geometry(crds_new_sf)[1:size_new], bg="red", cex=1.2 ,pch=21, add=TRUE) 	# out-of-sample points
  legend("bottomleft", legend=c("Data used for model","Data for prediction"), pch=c(21,21), pt.bg=c("blue","red"), bty="n")

  list(forecast.result = forecasts_result,
       RAMSE.pred = RAMSE_pred)
}

