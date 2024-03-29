#######################################
### Description of firms_sf dataset ###
#######################################
#
#' A collection of data on companies in the Lubelskie Voivodship in Poland.
#'
#' A collection of data on enterprises in the Lubelskie Voivodeship in Poland, based on the REGON register provided by the Polish Statistical Office
#' (Statistics Poland - Główny Urząd Statystyczny, GUS). It consists of a 10% representative sample of enterprises registered in 2012 in the NUTS2 region
#' of Lubelskie.
#'
#' @format An object of class `sf` with projection EPSG:4326 (WGS84) consisting of the following variables:
#' \item{empl}{Employment.}
#' \item{roa}{Return on assets.}
#' \item{dummy.agri}{Dummy variable indicating whether the enterprise belongs to the agricultural sector.}
#' \item{dummy.prod}{Dummy variable indicating whether the enterprise belongs to the production sector.}
#' \item{dummy.constr}{Dummy variable indicating whether the enterprise belongs to the construction sector.}
#' \item{dummy.serv}{Dummy variable indicating whether the enterprise belongs to the services sector.}
#' \item{dist.big.city}{Distance (in km) to the nearest large city.}
#' \item{dist.closest.city}{Distance (in km) to the nearest city.}
#' \item{year}{Year.}
#' \item{sector.aggr}{Sector to which the enterprise belongs.}
#' \item{sector}{Sector according to Polish Classification of Economic Activities (PKD 2007).}
#'
#' @source Own work based on REGON registry (Statistics Poland),
#' <https://wyszukiwarkaregon.stat.gov.pl/>
#'
#' @examples #To be done!!!
"firms_sf"
#
########################################
### Description of region_sf dataset ###
########################################
#
#' Boundaries of the Lubelskie Voivodship in Poland.
#'
#' Boundaries of the Lubelskie Voivodship in Poland (NUTS2 region) in object class `sf` obtained on the basis of shapefiles
#' of the Polish Boundary Register (PRG) provided by Head Office of Geodesy and Cartography in Poland
#' (Głowny Urząd Geodezji i Kartografii, GUGIK)
#'
#' @format An object of class `sf` with projection EPSG:4326 (WGS84) consisting of the boundaries of the Lubelskie Voivodeship in Poland.
#'
#' @source Based on data provided by Head Office of Geodesy and Cartography (GUIK),
#' <https://www.geoportal.gov.pl/pl/dane/panstwowy-rejestr-granic-prg/>
#'
#' @examples #To be done!!!
"region_sf"
