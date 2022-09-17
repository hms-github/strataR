#' Create a eustasy object
#'
#' Creates a eustasy object, which is needed to create a basin object. A eustasy object describes the history of eustatic sea level change through a basin simulation.
#' 
#' @title eustasy: The eustasy function
#' @param geometry a geometry object
#' @param period the period of eustatic sea-level cycles, in m.y.
#' @param amplitude the amplitude of eustatic sea-level cycles, in m. The peak-to-peak sea-level change is twice the amplitude.
#' @param symmetry a value ranging from 0 to 1 describing the symmetry of the eustatic cycles. 0.5 is a classic cosine wave; smaller values have increasingly faster falls and slower rises, and larger values have the opposite, slower falls and faster rises in sea level.
#' @param phase a string with possible values of rising, falling, highPoint, and lowPoint describing the starting position of eustatic sea level.
#' @param shape a value of 0 or greater that describes the shape of the eustatic cycle. 0 corresponds to a classic cosine curve, and larger values cause the curve to become increasingly squared-off.
#' @param netRise the net rise in meters of eustatic sea-level over the course of a basin simulation.
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, marginWidth=500, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' eust <- eustasy(geometry=geom, netRise=30.0)	
#' 
#' @rdname eustasy
#' @export eustasy

eustasy <- function(geometry, period=1, amplitude=0, symmetry=0.5, phase=c('rising', 'falling', 'highPoint', 'lowPoint'), shape=0, netRise=0) {
	# geometry is the object supplied by setGeometry()
	# all times (period, duration, timeStep, timePoint) are in m.y.
	# all lengths (amplitude, seaLevel, netRise) are in m
	# phase is one of four values - "falling", "rising", "highPoint", "lowPoint" - that describe whether the sine wave starts on the falling inflection point, rising inflection point, highest point, or lowest point on the sine wave.
	# shape is dimensionless, 0 is a sine wave, higher values are flattened sine waves
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

#' @return \code{NULL}
#' 
#' @rdname eustasy
#' @export

plot.eustasy <- function(x) {
	plot(x$timeSeries$timePoint, x$timeSeries$seaLevel, type='l', las=1, xlab="model time (m.y.)", ylab="sea level (m)")
}

#' @return \code{NULL}
#' 
#' @rdname eustasy
#' @export

print.eustasy <- function(x) {
	cat("net eustatic rise (netRise):                 ", x$parameters$netRise, "m\n")
	cat("Period of eustatic cyclicity (period):       ", x$parameters$period, "m.y.\n")
	cat("Amplitude of eustatic cyclicity (amplitude): ", x$parameters$amplitude, "m\n")
	cat("Symmetry of eustatic cyclicity (symmetry):   ", x$parameters$symmetry, "(dimensionless, 0 to 1)\n")
	cat("Phase of eustatic cyclicity (phase):         ", x$parameters$phase, "\n")
	cat("Shape of eustatic cyclicity (shape):         ", x$parameters$shape, "(dimensionless, 0 to infinity)\n")
	cat("Number of time steps:                        ", length(x$timeSeries$timePoint), "\n")
	cat("Minimum sea level:                           ", min(x$timeSeries$seaLevel), "m\n")
	cat("Maximum sea level:                           ", max(x$timeSeries$seaLevel), "m\n")
}

#' @return \code{NULL}
#' 
#' @rdname eustasy
#' @export

summary.eustasy <- function(x) {
	cat("net rise (netRise): ", x$parameters$netRise, "m\n")
	cat("period:             ", x$parameters$period, "m.y.\n")
	cat("amplitude:          ", x$parameters$amplitude, "m\n")
	cat("symmetry:           ", x$parameters$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:              ", x$parameters$phase, "\n")
	cat("shape:              ", x$parameters$shape, "(dimensionless, 0 to infinity)\n")
}



