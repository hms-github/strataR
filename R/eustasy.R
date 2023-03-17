#' @title Eustatic Sea-level History
#'
#' @description Create an object of class "eustasy".
#'
#' @details Creates a eustasy object, which is needed to create a basin object. A eustasy object describes the history of eustatic sea level change through a basin simulation.
#' 
#' @param geometry an object of class [`geometry`].
#' @param period the period of eustatic sea-level cycles, in m.y.
#' @param amplitude the amplitude of eustatic sea-level cycles, in m. The peak-to-peak  
#'   sea-level change is twice the amplitude.
#' @param symmetry a value ranging from 0 to 1 describing the symmetry of the eustatic  
#'   cycles. 0.5 is a classic cosine wave; smaller values have increasingly faster falls  
#'   and slower rises, and larger values have the opposite, slower falls and faster rises  
#'   in sea level.
#' @param phase a string with possible values of "falling", "rising", "highPoint", and  
#'   "lowPoint" describing the starting position of eustatic sea level.
#' @param shape a dimensionless value of 0 or greater that describes the shape of the  
#'   eustatic cycle. 0 corresponds to a classic cosine curve, and larger values cause the  
#'   curve to become increasingly squared-off.
#' @param netRise the net rise in meters of eustatic sea-level over the course of a basin  
#'   simulation.
#' @param x,object an object of class `eustasy`.
#' @param ... additional arguments to be passed.
#' 
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, 
#'   marginWidth=600, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' eust <- eustasy(geometry=geom, netRise=30.0)	
#' summary(eust)
#' print(eust)
#' plot(eust)
#' 
#' @rdname eustasy
#' @export eustasy
#' 
#' @return `eustasy` returns an object of class "eustasy", which includes print, summary, and plot methods.
#' 
#' A eustasy object consists of a two-item list. The first item is a list of the arugments that were used to create the eustatic history. The second item is a data frame called `timeSeries`, in which the first column is the model time point (`timePoint`) and the second is the sea level in meters (`seaLevel`).
#' 

eustasy <- function(geometry, period=1, amplitude=0, symmetry=0.5, phase=c('rising', 'falling', 'highPoint', 'lowPoint'), shape=0, netRise=0) {
	phase <- match.arg(phase)
	timePoint <- seq(0, geometry$duration, geometry$timeStep)
	linearTrend <- netRise * timePoint/geometry$duration
	seaLevel <- flexSin(timePoint, period=period, amplitude=amplitude, symmetry=symmetry, phase=phase, shape=shape) + linearTrend
	timeSeries <- data.frame(cbind(timePoint=timePoint, seaLevel=seaLevel))
	parameters <- list(netRise=netRise, period=period, amplitude=amplitude, symmetry=symmetry, phase=phase, shape=shape)
	results <- list(parameters=parameters, timeSeries=timeSeries)
	class(results) <- "eustasy"
	return(results)
}

#' @rdname eustasy
#' @export

plot.eustasy <- function(x, ...) {
	plot(x$timeSeries$timePoint, x$timeSeries$seaLevel, type='l', las=1, xlab="model time (m.y.)", ylab="sea level (m)", ...)
}

#' 
#' @rdname eustasy
#' @export

print.eustasy <- function(x, ...) {
	x
}

#' 
#' @rdname eustasy
#' @export

summary.eustasy <- function(object, ...) {
	cat("net rise (netRise):   ", object$parameters$netRise, "m\n")
	cat("period:               ", object$parameters$period, "m.y.\n")
	cat("amplitude:            ", object$parameters$amplitude, "m\n")
	cat("symmetry:             ", object$parameters$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:                ", object$parameters$phase, "\n")
	cat("shape:                ", object$parameters$shape, "(dimensionless, 0 to infinity)\n")
	cat("Number of time steps: ", length(object$timeSeries$timePoint), "\n")
	cat("Minimum sea level:    ", min(object$timeSeries$seaLevel), "m\n")
	cat("Maximum sea level:    ", max(object$timeSeries$seaLevel), "m\n")
}



