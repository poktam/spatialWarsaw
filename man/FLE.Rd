% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fle.R
\name{FLE}
\alias{FLE}
\alias{fle}
\title{Focal Local Entropy for rasterised data to analyse the density of points}
\usage{
FLE(points_sf, region_sf, nrows.raster = 50, ncols.raster = 50, w, r)

fle(points_sf, region_sf, nrows.raster = 50, ncols.raster = 50, w, r)
}
\arguments{
\item{points_sf}{Object in \code{sf} or the \code{data.frame} class - in the case of a \code{data.frame} object, the first and second columns must contain X and Y coordinates.
For \code{data.frame}, make sure that the coordinates of the points are in the same coordinate system / projection as the \code{region_sf} object.}

\item{region_sf}{The boundary of the area where the \code{points_sf} are located, \code{sf} class object.}

\item{nrows.raster}{Number of rows in the raster, default 50.}

\item{ncols.raster}{Number of columns in the raster, default 50.}

\item{w}{The size of the square moving window (choose an odd value such as 3,5,7,9 etc.); the w and r parameters are mutually exclusive - specify either w or r.}

\item{r}{The radius of a circular ring in a scale of geo-coordinates; the w and r parameters are mutually exclusive - specify either w or r.}
}
\value{
\code{FLE()} returns the \code{terra} class object and its visualisation.
}
\description{
The function analyses the density of geolocated points by counting the number of points in raster cells
and calculating local entropy using neighbouring cells. The counts in the grid cells are normalised
and the entropy is calculated for the discretised variable using fixed intervals. The focal mechanism performs
the calculations for each grid cell and its neighbours. The nearest neighbours can be given as a number
(an odd value e.g. 3,5,7,9) or a radius.
}
\details{
The function grids the area and counts the number of observations in each grid cell. The counts of the observations are normalised.
The focal mechanism selects the neighbourhood for each cell and creates a vector of normalised counts for the analysed cell
and its neighbours. The values are discretised into fixed intervals \code{(breaks=c(-100, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 100))}
and the proportion (percentage of observations) in each interval is calculated. These proportions are used as input to the entropy function.
The local entropy is calculated for each grid cell and displayed in the raster.

The choice of w or r makes a small difference to the result. For w=9, a cell and all its direct neighbours (first row, contiguity)
are considered. When setting r, remember that the distance to the centroids of neighbouring cells is closer
for east/west/north/south cells and longer for east-north, east-south, west-north and west-south cells.
Too short r can exclude the second group from the neighbours list.
}
\examples{
my.fle<-FLE(firms_sf, region_sf, nrows.raster=50, ncols.raster=50, w=9)
my.fle

}
