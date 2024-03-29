% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eustasy.R
\name{eustasy}
\alias{eustasy}
\alias{plot.eustasy}
\alias{summary.eustasy}
\title{Eustatic Sea-level History}
\usage{
eustasy(
  geometry,
  period = 1,
  amplitude = 0,
  symmetry = 0.5,
  phase = c("rising", "falling", "highPoint", "lowPoint"),
  shape = 0,
  netRise = 0
)

\method{plot}{eustasy}(x, ...)

\method{summary}{eustasy}(object, ...)
}
\arguments{
\item{geometry}{an object of class \code{\link{geometry}}.}

\item{period}{the period of eustatic sea-level cycles, in m.y.}

\item{amplitude}{the amplitude of eustatic sea-level cycles, in m. The peak-to-peak
sea-level change is twice the amplitude.}

\item{symmetry}{a value ranging from 0 to 1 describing the symmetry of the eustatic
cycles. 0.5 is a classic cosine wave; smaller values have increasingly faster falls
and slower rises, and larger values have the opposite, slower falls and faster rises
in sea level.}

\item{phase}{a string with possible values of "falling", "rising", "highPoint", and
"lowPoint" describing the starting position of eustatic sea level.}

\item{shape}{a dimensionless value of 0 or greater that describes the shape of the
eustatic cycle. 0 corresponds to a classic cosine curve, and larger values cause the
curve to become increasingly squared-off.}

\item{netRise}{the net rise in meters of eustatic sea-level over the course of a basin
simulation.}

\item{x, object}{an object of class \code{eustasy}.}

\item{...}{additional arguments to be passed.}
}
\value{
\code{eustasy} returns an object of class "eustasy", which includes print, summary, and plot methods.

A eustasy object consists of a two-item list. The first item is a list of the arugments that were used to create the eustatic history. The second item is a data frame called \code{timeSeries}, in which the first column is the model time point (\code{timePoint}) and the second is the sea level in meters (\code{seaLevel}).
}
\description{
Create an object of class "eustasy".
}
\details{
Creates a eustasy object, which is needed to create a basin object. A eustasy object describes the history of eustatic sea level change through a basin simulation.
}
\examples{
geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, 
  marginWidth=600, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
eust <- eustasy(geometry=geom, netRise=30.0)	
summary(eust)
plot(eust)

}
