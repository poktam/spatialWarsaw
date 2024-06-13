###########################
### ssr i pochodne    #####
###########################
#
#' @title Spatial Switching Regimes (SSR): estimation and comparison of econometric models in groups of different densities
#'
#' @description
#' The `ssr()` function works on geolocated point data. It uses input from the [spatialWarsaw::ClustConti()] or [spatialWarsaw::ClustDisjoint()] functions from the `spatialWarsaw` package to divide points
#' into density groups - within these subsamples it estimates the econometric model given by the equation. An important point is that the size and order of the datasets used to obtain density groups and
#' the econometric model must be the same.
#'
#' @details
#' Function returns a table with a summary of the models estimated in density groups.
#'
#' @name ssr
#' @param data_sf Geo-located points in `sf` or the `data.frame` class - in the case of a `data.frame` object, the first and second columns must contain X and Y coordinates.
#' ***NOTE! The same `data_sf` object must also be used to obtain a cluster split using the [spatialWarsaw::ClustConti()] or [spatialWarsaw::ClustDisjoint()] function.***
#' @param clusters Object of class `clust` as output of the function [spatialWarsaw::ClustConti()] or [spatialWarsaw::ClustDisjoint()] containing a vector of assignments of observations from the object `data_sf`
#' to individual clusters and the name of the function used for the division.
#' @param eq An object that defines the equation for the model, can be in the [stats::formula()] class (or one that can be coerced to that class) or a symbolic description of the model to be used.
#' @param family A `family` object  consistent with [stats::family()]: a description of the error distribution
#' and link function to be used in the model. The default value is `binomial` with logit link function.
#' @param ... Other parameters passed to [stats::glm()] during the calculation.
#'
#' @return `ssr()` returns the table with a summary of models estimated in density groups and a list object:
#' \item{data}{provides data by density groups.}
#' \item{coords}{provides all xy coordinates of point data, without the density group subdivision information.}
#' \item{cluster}{Provides a clustering vector that assigns all point data to density clusters.}
#' \item{type}{provides the type of clustering based on the `clusters` argument, one of "ClustDisjoint" or "ClustConti".}
#' \item{models}{privides the basic models (without significance tests) in each density cluster.}
#' The output also includes a `texreg` table which summarises the models in density groups and the size of density clusters.
#'
#' @examples
#' # The following example is on selected data from the firms_sf dataset included in the package
#' sub_sf<-firms_sf[1:5000,]
#' # One to choose
#' # clust<-ClustConti(sub_sf, clusters=4, noise =0.2, minPts0=30)
#' clust<-ClustDisjoint(sub_sf, noise=c(0.8, 0.6, 0.4, 0.2), minPts=30)
#'
#' # Density cluster plot
#' sub_sf$cluster<-clust$cluster
#' plot(sub_sf["cluster"], key.pos=4, pch=".", cex=2, main="Density clusters")
#'
#' # SSR models for the binary dependent variable
#' # The form of the equation to be estimated:
#' eq<-dummy.prod~empl+roa+dist.big.city
#'
#' my.ssr<-ssr(sub_sf, clust, eq, family=binomial)
#' my.ssr
#'
#' # Output with plot using ssrShowModels() function.
#' ssrShowModels(my.ssr, plot=TRUE)
#'
#' @seealso [ssrShowModels()]
#'
#' @export
ssr<-function(data_sf, clusters, eq, family=binomial, ...){

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

  if (!inherits(clusters,"clust")) {
    stop("The class of `clusters` argument must only be `clust`: result of ClustConti() or ClustDisjoint() from the spatialWarsaw package.")
  }

  if (nrow(data_sf)!=length(clusters$cluster)) {
    stop("The number of observations in the `data_sf` object and the length of the vector indicating
         cluster membership in the `clusters` object are not the same.")
  }

  var_names<-all.vars(eq)
  m <- match(gsub(" ", ".", var_names), colnames(data_sf))
  if (any(is.na(m))) {
    stop("Variable names in eq are incorrect.")
  }

  # czy nie sprawdzać jakoś parametru family -> to sprawdza glm

  data_sf_s<-split(data_sf, clusters$cluster) # create subgroups

  result_list<-list()
  for(i in 1:length(data_sf_s)){ # models in subgroups
    result_list[[i]]<-glm(formula=eq, family=family, data=data_sf_s[[i]], ...)
  }

  structure(
    list(
      data = data_sf_s, #czy to dawać? To są dane w podziale na subgrupy
      coords = crds,
      cluster = clusters$cluster,
      type = clusters$type,
      models = result_list
    ),
    class = "ssr"
  )

}

# The ssrShowModels() function is likely to be merged with the ssr() function in the future!!!!
#
#' @title Friendly summary and plot for ssr objects (ssr() function result)
#'
#' @description
#' The `ssrShowModels()` function with the `plot=TRUE` option shows the summary of the models and plots the points divided
#' by the applied density clusters.
#'
#' @details
#' The plot illustrates the division of geolocated points into density clusters. This partitioning is used in
#' the [ssr()] function for subgroup model estimation.
#'
#' @name ssrShowModels
#' @param ssr_model Result object of the [ssr()] function from this package. The class of `ssr_model` object must only be 'ssr'.
#' @param plot Logical; indicates whether the function should generate a plot, default `plot=TRUE`.
#' @param ... Other parameters passed to [stargazer::stargazer()].
#'
#' @return This function returns the summary table of econometric models in density clusters. If `plot=TRUE` it will also return
#' the the point plot with division into clusters.
#'
#' @examples
#' # The following example is on selected data from the firms_sf dataset included in the package
#' sub_sf<-firms_sf[1:5000,]
#' # One to choose
#' # clust<-ClustConti(sub_sf, clusters=4, noise =0.2, minPts0=30)
#' clust<-ClustDisjoint(sub_sf, noise=c(0.8, 0.6, 0.4, 0.2), minPts=30)
#'
#' # Density cluster plot
#' sub_sf$cluster<-clust$cluster
#' plot(sub_sf["cluster"], key.pos=4, pch=".", cex=2, main="Density clusters")
#'
#' # SSR models for the binary dependent variable
#' # The form of the equation to be estimated:
#' eq<-dummy.prod~empl+roa+dist.big.city
#'
#' my.ssr<-ssr(sub_sf, clust, eq, family=binomial)
#' my.ssr
#'
#' # Output with plot using ssrShowModels() function.
#' ssrShowModels(my.ssr, plot=TRUE)
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
    if(ssr_model$type == 'ClustConti'){
      cols = ifelse(ssr_model$cluster==0, max(ssr_model$cluster)+1, ssr_model$cluster)
      ggplot2::ggplot()+
        ggplot2::geom_point(ggplot2::aes(x=as.numeric(ssr_model$coords[,1]), y=as.numeric(ssr_model$coords[,2]), colour=as.factor(cols)))+
        xlab("")+
        ylab("")+
        ggplot2::ggtitle("Group assignments")+
        ggplot2::scale_color_viridis_d(labels = c(paste0("Cluster ", 1:max(ssr_model$cluster)),"Noise points"))+
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
  cat("\nAvailable fields: data, coords, cluster, type, models\n")
}
