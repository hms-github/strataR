% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsidence.R
\name{subsidence}
\alias{subsidence}
\alias{plot.subsidence}
\alias{summary.subsidence}
\title{Subsidence History}
\usage{
subsidence(
  geometry,
  startingLeft = 0,
  startingRight,
  netChangeFactor = 1,
  period = 1,
  amplitude = 0,
  symmetry = 0.5,
  phase = c("rising", "falling", "highPoint", "lowPoint"),
  shape = 0
)

\method{plot}{subsidence}(x, geometry, type = c("lines", "filled"), ...)

\method{summary}{subsidence}(object, ...)
}
\arguments{
\item{geometry}{an object of class \code{\link{geometry}}.}

\item{startingLeft}{the starting subsidence rate at the left edge of the model, where
sediment is introduced, in m/m.y.}

\item{startingRight}{the starting subsidence rate at the right edge of the model,
farthest from where sediment is introduced, in m/m.y.}

\item{netChangeFactor}{the factor by which subsidence changes from the beginning of the
simulation to the end. For example, 0.5 would indicate a halving of subsidence rates,
1 would indicate no change in rates, and 2 would indicate a doubling.}

\item{period}{the period of subsidence cycles, in m.y.}

\item{amplitude}{the amplitude of subsidence cycles, in m. The peak-to-peak subsidence
change is twice the amplitude.}

\item{symmetry}{a value ranging from 0 to 1 describing the symmetry of the subsidence
cycles. 0.5 is a classic cosine wave; smaller values have increasingly faster initial
change in rates and slower subsequent changes in rates; larger values indicate the
opposite.}

\item{phase}{a string with possible values of "rising", "falling", "highPoint", and
"lowPoint" describing the starting position of cyclic subsidence.}

\item{shape}{a dimensionless value of 0 or greater that describes the shape of the
subsidence cycle. 0 corresponds to a classic cosine curve, and larger values cause
the curve to become increasingly squared-off.}

\item{x, object}{an object of class \code{subsidence}.}

\item{type}{a string ("lines" or "filled") specifying whether contours should be shown
as lines or by a color gradient fill.}

\item{...}{additional arguments to be passed.}
}
\value{
subsidence returns an object of class "subsidence", which includes print, summary, and plot methods.

A subsidence object is a list of two items. The first item (\code{parameters}) is a list of the arguments used in creating the subsidence object. The second item (\code{rates}) is a matrix that records the subsidence rate (in m / m.y.) at every location in the basin at every time point. Rows correspond to time points, and columns correspond to positions in the basin.
}
\description{
Create an object of class "subsidence".
}
\details{
Creates a subsidence object, which is needed to create a basin object. A subsidence object describes how subsidence changes across the sedimentary basin and through time.
}
\examples{
geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, 
  marginWidth=600, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
subs <- subsidence(geometry=geom, startingLeft=0.0, startingRight=1.0)

}
