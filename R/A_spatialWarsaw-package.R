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
#' @importFrom benford.analysis benford
#' @importFrom dbscan dbscan kNN kNNdist frNN
#' @importFrom DescTools PseudoR2
#' @importFrom dplyr summarise group_by
#' @importFrom stats binomial glm dist rnorm kmeans
#' @importFrom sp coordinates
#' @importFrom sf st_as_sf st_bbox st_crs st_sample st_coordinates st_geometry_type st_voronoi st_cast st_union st_intersection st_area st_geometry st_transform
#' @importFrom stargazer stargazer
#' @importFrom ggplot2 ggplot geom_point aes aes_string ggtitle scale_color_viridis_d guides guide_legend theme_minimal theme element_text xlab ylab
#' @importFrom graphics legend par points
#'
# ew. uzupełnić inne importy!!!; może trzeba się pozbyć pakietu sp?
# SPRAWDZIĆ @importFrom vs Imports w DESCRIPTION i funkcjami użytymi
#'
NULL
