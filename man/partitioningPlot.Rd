% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partitioningPlot.R
\name{partitioningPlot}
\alias{partitioningPlot}
\title{Sediment Partitioning Plot}
\usage{
partitioningPlot(basin)
}
\arguments{
\item{basin}{an object of class \code{\link{basin}}.}
}
\description{
Plot showing how sediment is partitioned into nonmarine and marine areas through time.
}
\details{
Plots shows the volume of sediment deposited in nonmarine vs. marine areas through time, as a fraction of the total amount of sediment. As the shore transgresses and regresses, owing to subsidence, eustasy, and sediment supply, deposition may occur preferentially in nonmarine or marine areas.
}
\examples{
data(sedBasin)
partitioningPlot(sedBasin)

}
