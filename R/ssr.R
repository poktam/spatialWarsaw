###########################
### ssr i pochodne    #####
###########################
#
#' @title Linijka nr 1 - function title
#'
#' @description
#' Linijka nr 2 - description (razem do funkcji ssr i show_models?)
#'
#' @details
#' Linijka nr 3 - details (razem do funkcji ssr i show_models?)
#'
#' @name ssr
#' @param data_sf (obiekt 'sf' lub 'data.frame' zawierający współrzędne - w data frame 1 kolumna musi być X coords, druga kolumna Y coords)
#' @param type do opisu: one of: "ClustConti","ClustDisjoint".
#' @param clusters do opisu
#' @param sep do opisu
#' @param r_p do opisu
#' @param eps_r do opisu
#' @param eps_np do opisu
#' @param minPts do opisu
#' @param bootstrap do opisu
#' @param sample_size do opisu (used only if bootstrap=TRUE)
#' @param times do opisu  (used only if bootstrap=TRUE)
#' @param eq an object of class [stats::formula()] (or one that can be coerced to that class):
#' a symbolic description of the model to be used.
#' @param family a family object  consistent with [stats::family()]: a description of the error distribution
#' and link function to be used in the model. The default value is `binomial` with logit link function.
#' @param ... other parameters (do opisu - chyba przekazywane do glm())
#'
#' @return description
#'
#' @examples #To be done!!!
#'
#' @export
ssr<-function(data_sf, type="ClustConti", clusters, sep, r_p=0.001, eps_r=10e-16, eps_np=10e-3, minPts=5,
              bootstrap=FALSE, sample_size=NULL, times=NULL, eq, family=binomial, ...){

  params <- (as.list(environment()))[-1] # tu można poprawić, żeby w zależności od bootstrap dawało odpowiednie

  if((inherits(data_sf,"sf") && st_geometry_type(data_sf,FALSE)=="POINT")) {
    crds<-as.data.frame(st_coordinates(data_sf))
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class sf.\n", sep = "")
  }  else if(inherits(data_sf,"data.frame",TRUE)==1){
    crds<-data_sf[,c(1,2)]
    colnames(crds)<-c("X_coord","Y_coord")
    cat("Data_sf was detected as an object of class data.frame.\n", sep = "")
  }  else {
    stop("The class of data_sf must only be 'sf' of geometry type 'POINTS' or 'data.frame'.")
  }

  # zbadać type
  if (length(type)>1) {
    type<-type[1]
    cat("Parameter type longer than 1. The first element has been selected: ",type,"\n",sep="")
  }

  if (!(type %in% c("ClustConti","ClustDisjoint"))) {
    stop("Unknown model type. Must be one of: ClustConti, ClustDisjoint")
  }

  # zbadać clusters, czy nie ujemne - może inny warunek niż <=0
  if (type=="ClustConti"){
    clusters<-round(clusters,0)
    if (clusters<=0) {
      clusters<-4
      cat("Incorrect clusters parameter. Set to the proposed value: 4.","\n",sep="")
    }
  }

  # zbadać minPts, czy nie ujemne - może inny warunek niż <=0
  minPts<-round(minPts,0)
  if (minPts<=0) {
    minPts<-5
    cat("Incorrect minPts parameter. Set to the proposed value: 5.","\n",sep="")
  }

  if(bootstrap){

    if (is.null(sample_size) || is.null(times)) {
      stop("For bootstrap=TRUE, both sample_size and times parameters must be set.")
    }

    # zbadać sample_size i ustawić ew. na wielkość zbioru danych
    sample_size <- as.integer(sample_size)
    if (sample_size > nrow(crds) || sample_size<1) {
      sample_size<-nrow(crds)
      cat("Wrong sample size. Sample_size set to:",sample_size,"\n",sep="")
    }

    # zbadać times, czy nie ujemne - może inny warunek niż <=0
    times<-round(times,0)
    if (times<=0) {
      times<-16
      cat("Incorrect times value. Times set to the proposed value: 16.","\n",sep="")
    }
  }

  var_names<-all.vars(eq)
  m <- match(gsub(" ", ".", var_names), colnames(data_sf))
  if (any(is.na(m))) {
    stop("Variable names in eq are incorrect.")
  }

  # czy nie sprawdzać jakoś parametru family -> to sprawdza glm

  if(type=="ClustConti" & bootstrap==FALSE)
  { # group points
    output<-ClustConti(data_sf=crds, clusters=clusters, noise=sep, r_p=r_p, eps_r=eps_r, eps_np=eps_np)$cluster
  } else if(type=="ClustConti" & bootstrap==TRUE) {
    stop("Bootstraped DBSCAN cannot be applied to the ClustConti version of this model")
  } else if(type=="ClustDisjoint" & bootstrap==FALSE) {
    output<-ClustDisjoint(data_sf=crds, noise=sep, r_p=r_p, eps_r=eps_r, eps_np=eps_np, minPts=minPts,
                          bootstrap=FALSE, sample_size=NULL, times=NULL)$cluster
  } else if(type=="ClustDisjoint" & bootstrap==TRUE) {
    output<-ClustDisjoint(data_sf=crds, noise=sep, r_p=r_p, eps_r=eps_r, eps_np=eps_np, minPts=minPts,
                          bootstrap=TRUE, sample_size=sample_size, times=times)$cluster
  }

  data_sf_s<-split(data_sf, output) # create subgroups

  result_list<-list()
  for(i in 1:length(data_sf_s)){ # models in subgroups
    result_list[[i]]<-glm(formula=eq, family=family, data=data_sf_s[[i]], ...)
  }

  structure(
    list(
      data = data_sf_s, #czy to dawać? To są dane w podziale na subgrupy
      coords = crds,
      models = result_list,
      cluster = output,
      parameters = params
    ),
    class = "ssr"
  )

}

#' @title Linijka nr 1 - ssrShowModels() title
#'
#' @description
#' A short description of ssrShowModels
#'
#' @details
#' Linijka nr 3 - details (razem do funkcji ssr i show_models?)
#'
#' @name ssrShowModels
#' @param ssr_model do opisu
#' @param plot do opisu
#'
#' @return ssrShowModels description (DZIAŁA)
#'
#' @examples #To be done!
#'
#' @export
ssrShowModels<-function(ssr_model, plot=TRUE, ...){
  if (!inherits(ssr_model,"ssr")) {
    stop("The class of ssr_model must only be 'ssr': result of ssr() from spatialWarsaw package.")
  }

  psR2<-lapply(ssr_model$models, function(x) DescTools::PseudoR2(x, which="all"))
  R2mcf<-c("R2 McFadden", round(sapply(psR2, "[[", 1),3))
  R2n<-c("R2 Nagelkerke", round(sapply(psR2, "[[", 4),3))
  R2mkz<-c("R2 McKelvey.Zavoina", round(sapply(psR2, "[[", 8),3))
  labels<-c("Most dense", rep("", length(ssr_model$models)-2),"Least dense")
  stargazer(ssr_model$models, type="text", add.lines=list(R2mcf, R2n, R2mkz), column.labels = labels, ...)

  if(plot){
    if(ssr_model$parameters$type == 'ClustConti'){
      cols = ifelse(ssr_model$cluster==0, ssr_model$parameters$clusters+1, ssr_model$cluster)
      ggplot2::ggplot()+
        ggplot2::geom_point(ggplot2::aes(x=as.numeric(ssr_model$coords[,1]), y=as.numeric(ssr_model$coords[,2]), colour=as.factor(cols)))+
        xlab("")+
        ylab("")+
        ggplot2::ggtitle("Group assignments")+
        ggplot2::scale_color_viridis_d(labels = c(paste0("Cluster ", 1:ssr_model$parameters$clusters),"Noise points"))+
        ggplot2::guides(col=ggplot2::guide_legend(title="Group:"))+
        ggplot2::theme_minimal()+
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
      }
    else{
      ggplot2::ggplot()+
        ggplot2::geom_point(ggplot2::aes(x=as.numeric(ssr_model$coords[,1]), y=as.numeric(ssr_model$coords[,2]), colour=as.factor(ssr_model$cluster)))+
        ggplot2::xlab("")+
        ggplot2::ylab("")+
        ggplot2::ggtitle("Group assignments")+
        ggplot2::scale_color_viridis_d(labels = c("Most\ndense",seq(2, length(ssr_model$models)-1),"Least\ndense"))+
        ggplot2::guides(col=ggplot2::guide_legend(title="Group:"))+
        ggplot2::theme_minimal()+
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5))
      }
  }
}


#' @export
# define what is the output of print(ssr)
# MK: uproszczone względem oryginalnej wersji, można wrócić do starej jeśli nie działa.
print.ssr <- function(x){
#  cat("Call:\n",x$this.call,"\n")
  cat("Group assignments:\n",summary(factor(x$cluster)),"\n")
  cat("\nSummary of model coefficients:\n")
  cat(ssrShowModels(x, plot=FALSE))
  cat("\nAvailable fields: data, coords, models, cluster, parameters\n")
}
