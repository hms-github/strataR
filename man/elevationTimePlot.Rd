% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elevationTimePlot.R
\name{elevationTimePlot}
\alias{elevationTimePlot}
\title{Elevation–Time Plot}
\usage{
elevationTimePlot(basin, locationKm, setting = c("valley", "interfluve"), ...)
}
\arguments{
\item{basin}{an object of class \code{\link{basin}}.}

\item{locationKm}{the location in km as a distance from the left edge of the basin}

\item{setting}{a string, either "valley" or "interfluve".}

\item{...}{Optional arguments to be passed to the plot.}
}
\description{
Plot of elevation through time.
}
\details{
Plot showing how the elevation (in m above sea level) at a location in the basin changes over time.
}
\examples{
data(sedBasin)
elevationTimePlot(sedBasin, locationKm=100, setting="valley")

}
