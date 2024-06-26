% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootspatreg.R
\name{SpatPredTess}
\alias{SpatPredTess}
\title{Prediction for out-of-sample point with spatial econometric models using Voronoi tesselation}
\usage{
SpatPredTess(model_spatial, points_spatial_sf, knnW, points_new_sf, region_sf)
}
\arguments{
\item{model_spatial}{Spatial model object to predict with.}

\item{points_spatial_sf}{Geo-located points in \code{sf} or the \code{data.frame} class used in the estimation of \code{model_spatial} -
in the case of a \code{data.frame} object, the first and second columns must contain X and Y coordinates.
For \code{data.frame}, make sure that the coordinates of the points are in the same coordinate system / projection as the \code{region_sf} object.}

\item{knnW}{Spatial weighting matrix (\code{listw} class) used to estimate \code{model_spatial}.}

\item{points_new_sf}{New data for the forecast. The new data must have the same class, structure and projection as the \code{points_spatial_sf} object.}

\item{region_sf}{Polygon in the \code{sf} class that defines the boundary for \code{points_spatial_sf}. New points from \code{points_new_sf} must be inside this boundary.}
}
\value{
\code{SpatPredTess()} returns the characteristics of new prediction points: values of predicted and real y, x and y coordinates of the point, and squared distance between real and predicted values (predY-realY)^2.

\code{SpatPredTess()} also returns the Voronoi diagram of the in-sample points with red new points and the quality of the prediction RMSE.
}
\description{
Spatial econometric models estimated on point data are limited in predicting values for the new observations that are not included
in the spatial weight matrix W. This function allows predictions for out-of-sample point data beyond W, as long as they are within
the envelope of existing points. This algorithm performs a Voronoi tesselation of the in-sample points and checks to which
Voronoi tile the new out-of-sample point belongs. It then replaces the data of the original point, which was assigned
to the same Voronoi tile, with the values of a new point. Predictions for a new point are made using n-1 old points
(which help to build spatial lags) and one new point.
}
\details{
Spatial predictions are only made for a single new point. If the input specifies more than one new point,
the predictions are made sequentially, e.g. for 10 new points the algorithm will make 10 prediction runs.
}
\examples{
# sf and spdep packages required
library(sf)
library(spdep)

# The input to SpatPredTess() is the object with the estimated model.
# The first two lines are taken from the BootSpatReg() example.
eq<-roa~empl+dummy.prod+dummy.constr+dummy.serv+dist.big.city
bsr<-BootSpatReg(firms_sf, 5, 1500, eq, "SEM", knn=5)

# Coordinates in sf class - observations used only in the best model.
# Required for knn spatial weighting matrix.
crds.sf<-st_centroid(st_geometry(firms_sf[rownames(bsr$data.best),]))

# Spatial weight matrix for predictive model: symmetric k nearest neighbours (knn=6) matrix used
knn.sym.listw<-nb2listw(make.sym.nb(knn2nb(knearneigh(crds.sf, k=6))))

# Prediction using Voronoi tiles
pred<-SpatPredTess(bsr$model.best, firms_sf[rownames(bsr$data.best),], knn.sym.listw, firms_sf[1:10,], region_sf)
pred


}
\references{
Kopczewska, K. (2023). Spatial bootstrapped microeconometrics: Forecasting for out‐of‐sample geo‐locations in big data. Scandinavian Journal of Statistics.
}
