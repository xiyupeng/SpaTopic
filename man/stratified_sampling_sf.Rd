% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_func.R
\name{stratified_sampling_sf}
\alias{stratified_sampling_sf}
\title{Spatially stratified random sample points from an image.}
\usage{
stratified_sampling_sf(
  points,
  cellsize = c(600, 600),
  num_samples_per_stratum = 1
)
}
\arguments{
\item{points}{a data frame contains all points in a image with X, Y coordinates.}

\item{cellsize}{a vector of length 2 contains the size of each grid square. 
Default c(600,600).}

\item{num_samples_per_stratum}{number of point selected from each grid square.
Default 1.}
}
\value{
Return a vector contains index of sampled points.
}
\description{
Spatially stratified random sample points from an image by R package \code{sf}
}
\examples{

data("lung5")
pt_idx<-stratified_sampling_sf(lung5, cellsize = c(600,600))

}
