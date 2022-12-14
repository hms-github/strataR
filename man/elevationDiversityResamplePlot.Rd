% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elevationDiversityResamplePlot.R
\name{elevationDiversityResamplePlot}
\alias{elevationDiversityResamplePlot}
\title{Resampled Elevation-Diversity Plot}
\usage{
elevationDiversityResamplePlot(
  species,
  minElevation = 0,
  maxElevation = 150,
  time = 0,
  numTrials = 100,
  conf = 0.5
)
}
\arguments{
\item{species}{an object of class \code{\link{species}}.}

\item{minElevation, maxElevation}{a numeric value for the minimum and maximum elevation
to be considered, in meters above sea level.}

\item{time}{a numeric value, in m.y., that limits analysis to species that originated
before a specified time and that went extinct after that same time. Generally left at
default to include species alive at beginning of simulation.}

\item{numTrials}{a numeric value for the number of resampling trials.}

\item{conf}{is a numeric value (0 to 1) indicating the confidence level used for
drawing the uncertainty envelope}
}
\description{
A plot of diversity vs. elevation based on resampling.
}
\details{
Diagnostic plot to determine if diversity covaries with elevation; it shouldn't, but this might happen for stochastic reasons, especially if the number of species is small. Plot spans a specified range of elevations and resamples to test if each species is present at a particular elevation. See \code{elevationDiversityPlot()} for a faster version that lacks confidence intervals and is based on a fixed threshold probability of collection.
}
\examples{
data(spec)
elevationDiversityResamplePlot(species=spec, minElevation=0, maxElevation=150, 
  numTrials=100) 

}
