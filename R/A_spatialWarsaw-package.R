#' @title `r packageDescription("spatialWarsaw")$Package`: `r packageDescription("spatialWarsaw")$Title`
#'
#' @description `r sub("<[^>]+>", "", packageDescription("spatialWarsaw")$Description)`
#'
#' @author Katarzyna Kopczewska and Mateusz Kopyt
#'
#' @references
#' To be done
#'
#' @name spatialWarsaw-package
#'
#' @aliases spatialWarsaw-package NULL
#' @keywords internal
"_PACKAGE"
#'
#' @import sf
#' @importFrom benford.analysis benford
#' @importFrom dbscan dbscan kNN kNNdist frNN
#' @importFrom DescTools PseudoR2
#' @importFrom dplyr summarise group_by
#' @importFrom stats binomial glm dist rnorm kmeans median AIC sd cor uniroot
#' @importFrom sp coordinates
#' @importFrom stargazer stargazer
#' @importFrom ggplot2 ggplot geom_point aes aes_string ggtitle scale_color_viridis_d guides guide_legend theme_minimal theme element_text xlab ylab
#' @importFrom graphics legend par points abline text title
#' @importFrom rlang check_exclusive
#' @importFrom terra rast rasterize focalMat focal values plot scale vect buffer aggregate expanse
#' @importFrom GWmodel bw.gwr gwr.basic gw.dist
#' @importFrom fossil adj.rand.index
#' @importFrom grDevices heat.colors
#' @importFrom spdep poly2nb nb2listw knearneigh make.sym.nb knn2nb lag.listw nb2mat
#' @importFrom spatialreg lagsarlm errorsarlm sacsarlm
#' @importFrom gridExtra grid.arrange
#'
#'
# ew. uzupełnić inne importy!!!; sp trzeba na razie zostawić, bo GWR models działają na sp
# może dla plot.matrix spróbować: @importMethodsFrom package generic !!!! - pozbycie się Depend z DESCRIPTION
#
#'
NULL
