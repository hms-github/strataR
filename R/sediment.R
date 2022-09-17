#' Create a sediment object
#'
#' Creates a sediment object, which is needed to create a basin object. A sediment object describes how sediment flux to the basin changes through time
#' 
#' @title sediment: The sediment function
#' @param geometry a geometry object
#' @param startingVolume the starting sediment flux delivered at the left edge of the model, in m * km.
#' @param netIncrease the amount by which sediment influx increases beginning of the simulation to the end. Negative values indicate a decrease in sediment flux.
#' @param period the period of sediment flux cycles, in m.y.
#' @param amplitude the amplitude of sediment flux cycles, in m * km. The peak-to-peak sediment flux change is twice the amplitude.
#' @param symmetry a value ranging from 0 to 1 describing the symmetry of the sediment flux cycles. 0.5 is a classic cosine wave; smaller values have increasingly faster initial change in rates and slower subsequent changes in rates; larger values indicate the opposite.
#' @param phase a string with possible values of rising, falling, highPoint, and lowPoint describing the starting position of cyclic sediment flux.
#' @param shape a value of 0 or greater that describes the shape of the sediment flux cycle. 0 corresponds to a classic cosine curve, and larger values cause the curve to become increasingly squared-off.
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, marginWidth=500, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' sedi <- sediment(geometry=geom, startingVolume=60)
#' 
#' @rdname sediment
#' @export sediment

sediment <- function(geometry, startingVolume=60, netIncrease=0, period=1, amplitude=0, symmetry=0.5, phase=c('rising', 'falling', 'highPoint', 'lowPoint'), shape=0) {
	# geometry is the object supplied by setGeometry()
	# all times (period, duration, timeStep, timePoint) are in m.y.
	# startingVolume, netIncrease, and amplitude are sediment volumes
	# phase is one of four values - "falling", "rising", "highPoint", "lowPoint" - that describe whether the sine wave starts on the falling inflection point, rising inflection point, highest point, or lowest point on the sine wave.
	# shape is dimensionless, 0 is a sine wave, higher values are flattened sine waves
	phase <- match.arg(phase)
	timePoint <- seq(0, geometry$duration, geometry$timeStep)
	linearTrend <- netIncrease * timePoint/geometry$duration
	sedimentVolume <- startingVolume + linearTrend + flexSin(timePoint, period=period, amplitude=amplitude, symmetry=symmetry, phase=phase, shape=shape)
	if (min(sedimentVolume) <= 0) {
		warning('zero or negative sedimentation rates generated, replacing with a zero; should redesign sediment generation to avoid this', call.=FALSE, immediate.=TRUE, noBreaks.=TRUE) 
		sedimentVolume[sedimentVolume <= 0] <- 0.0
	}
	timeSeries <- data.frame(cbind(timePoint=timePoint, volume=sedimentVolume))
	parameters <- list(startingVolume=startingVolume, netIncrease=netIncrease, period=period, amplitude=amplitude, symmetry=symmetry, phase=phase, shape=shape)
	results <- list(parameters=parameters, timeSeries=timeSeries)
	class(results) <- "sediment"
	return(results)
}

#' @return \code{NULL}
#' 
#' @rdname sediment
#' @export

plot.sediment <- function(x) {
	plot(x$timeSeries$timePoint, x$timeSeries$volume, type='l', las=1, xlab="model time (m.y.)", ylab="sediment flux")
}

#' @return \code{NULL}
#' 
#' @rdname sediment
#' @export

print.sediment <- function(x) {
	cat("startingVolume:          ", x$parameters$startingVolume, "m^2*km\n")
	cat("netIncrease:             ", x$parameters$netIncrease, "m^2*km\n")
	cat("period:                  ", x$parameters$period, "m.y.\n")
	cat("amplitude:               ", x$parameters$amplitude, "m^2*km\n")
	cat("symmetry:                ", x$parameters$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:                   ", x$parameters$phase, "\n")
	cat("shape:                   ", x$parameters$shape, "(dimensionless, 0 to infinity)\n")
	cat("timePoints:              ", length(x$timeSeries$timePoint), "\n")
	cat("minimum sedment volume:  ", min(x$timeSeries$volume), "m^2*km\n")
	cat("maximum sediment volume: ", max(x$timeSeries$volume), "m^2*km\n")
}

#' @return \code{NULL}
#' 
#' @rdname sediment
#' @export

summary.sediment <- function(x) {
	cat("Starting volume of sediment (startingVolume):     ", x$parameters$startingVolume, "m^2*km\n")
	cat("Net increase in sediment (netIncrease):           ", x$parameters$netIncrease, "m^2*km\n")
	cat("Period of cyclical sediment input (period):       ", x$parameters$period, "m.y.\n")
	cat("Amplitude of cyclical sediment input (amplitude): ", x$parameters$amplitude, "m^2*km\n")
	cat("Symmetry of cyclical sediment input (symmetry):   ", x$parameters$symmetry, "(dimensionless, 0 to 1)\n")
	cat("Phase of cyclical sediment input (phase):         ", x$parameters$phase, "\n")
	cat("Shape of cyclical sediment input (shape):         ", x$parameters$shape, "(dimensionless, 0 to infinity)\n")
}
