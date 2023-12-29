#' @title `r packageDescription("spatialWarsaw")$Package`: `r packageDescription("spatialWarsaw")$Title`
#'
#' @description `r sub("<[^>]+>", "", packageDescription("spatialWarsaw")$Description)`
#'
#' @author Katarzyna Kopczewska and Mateusz Kopyt
#'
#' @references
#' To be done
#'
#' @docType package
#' @name spatialWarsaw-package
#'
#' @import sf
#' @importFrom benford.analysis benford
#' @importFrom dbscan dbscan kNN kNNdist frNN
#' @importFrom DescTools PseudoR2
#' @importFrom dplyr summarise group_by
#' @importFrom stats binomial glm dist rnorm kmeans
#' @importFrom sp coordinates
#' @importFrom stargazer stargazer
#' @importFrom ggplot2 ggplot geom_point aes aes_string ggtitle scale_color_viridis_d guides guide_legend theme_minimal theme element_text xlab ylab
#' @importFrom graphics legend par points
#' @importFrom rlang check_exclusive
#' @importFrom terra rast rasterize focalMat focal values plot
#'
# ew. uzupełnić inne importy!!!; może trzeba się pozbyć pakietu sp?
# SPRAWDZIĆ @importFrom vs Imports w DESCRIPTION i funkcjami użytymi
#'
NULL
