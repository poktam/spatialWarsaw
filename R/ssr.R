###########################
### ssr i pochodne    #####
###########################
### DO TOTALNEJ POPRAWY ###
###########################
# MK: ustalić czy ssr ma być dostępna zewnątrz, ale chyba tak

#' Linijka nr 1 - function title
#'
#' Linijka nr 2 - description (razem do funkcji ssr i show_models)
#'
#' Linijka nr 3 - details (razem do funkcji ssr i show_models)
#'
#'
#' @name ssr
#' @param formula do opisu
#' @param data.spsf do opisu
#' @param type do opisu
#' @param clusters do opisu
#' @param minPts do opisu
#' @param r_p do opisu
#' @param sep do opisu
#' @param bootstrap do opisu
#' @param sample_size do opisu
#' @param times do opisu
#' @param family do opisu
#' @param eps_r do opisu
#' @param eps_np do opisu
#' @param ... other paramethers (do opisu)
#' @examples #To be done!!!
#' @export
# może parametr family uprościć, żeby podawać tylko nazwę funkcji łączącej
ssr<-function(formula, data.spsf, type, clusters, minPts=5, r_p=0.001, sep, bootstrap=FALSE, sample_size, times,
              family=stats::binomial(link="logit"), eps_r=10e-16, eps_np=10e-3, ...){

  data.df<-as.data.frame(data.spsf)

  if("sf" %in% class(data.spsf)) {crds=sp::coordinates(data.spsf)}
  else if("sf" %in% class(data.spsf)){crds=sf::st_coordinates(data.spsf)}
  else {
    stop("The class of data.spsf must be 'sp' or 'sf'.") # sp zachowane tymczasowo dla wstecznej zgodności
  }

  if(type=="A" & bootstrap==0) { # group points
    output<-ClustConti(crds, clusters, sep, r_p, eps_r, eps_np)$cluster }else if(type=="A" & bootstrap==1) {
      stop("Bootstraped DBSCAN cannot be applied to the A version of this model")}else if(type=="B" & bootstrap==0) {
        output<-ClustDisjoint(crds, minPts, sep, r_p, eps_r, eps_np, bootstrap=FALSE, 0, 0)$cluster}else if(type=="B" & bootstrap==1) {
          output<-ClustDisjoint(crds, minPts, sep, r_p, eps_r, eps_np, bootstrap=TRUE, sample_size, times)$cluster}

  data_s<-split(data.df, output) # create subgroups

  result_list<-list()
  for(i in 1:length(data_s)){ # models in subgroups
    result_list[[i]]<-stats::glm(formula=formula, data=data_s[[i]], family=family, ...)
  }

  structure(
    list(
      data = data.df,
      coords = crds,
      models = result_list,
      cluster = output,
      this.call = match.call()
    ),
    class = "ssr"
  )

}

# show_models MK: do sprawdzenia
# MK: czy dostępna z zewnątrz
# MK: z funkcją ggplot2::aes() może być problem - użycie aes w pakietach do sprawdzenia
#' A TU MOŻE TROSZKĘ POMOCY DO FUNKCJI show_models jeśli dostępne na zewnątrz :)
#' @rdname ssr
#' @param ssr_model do opisu
#' @param plot do opisu
#' @export
show_models<-function(ssr_model, plot=TRUE, ...){
  psR2<-lapply(ssr_model$models, function(x) DescTools::PseudoR2(x, which="all"))
  R2mcf<-c("R2 McFadden", round(sapply(psR2, "[[", 1),3))
  R2n<-c("R2 Nagelkerke", round(sapply(psR2, "[[", 4),3))
  R2mkz<-c("R2 McKelvey.Zavoina", round(sapply(psR2, "[[", 8),3))
  labels = c("Most dense", rep("", length(ssr_model$models)-2),"Least dense")
  stargazer::stargazer(ssr_model$models, type="text", add.lines=list(R2mcf, R2n, R2mkz), column.labels = labels, ...)

  if(plot){
    if(ssr_model$this.call$type == 'A'){
      cols = ifelse(ssr_model$cluster==0, ssr_model$this.call$clusters+1, ssr_model$cluster)
      ggplot2::ggplot()+
        ggplot2::geom_point(ggplot2::aes(x=as.numeric(ssr_model$coords[,1]), y=as.numeric(ssr_model$coords[,2]), colour=as.factor(cols)))+
        xlab("")+
        ylab("")+
        ggplot2::ggtitle("Group assignments")+
        ggplot2::scale_color_viridis_d(labels = c(paste0("Cluster ", 1:ssr_model$this.call$clusters),"Noise points"))+
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
  cat("Call:\n",x$this.call,"\n")
  cat("Group assignments:\n",summary(factor(x$cluster)),"\n")
  cat("Summary of model coefficients:\n")
  cat(show_models(x, plot=FALSE))
}
