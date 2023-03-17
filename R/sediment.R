#' @title Sediment Flux History
#' 
#' @description Create an object of class "sediment".
#'
#' @details Creates a sediment object, which is needed to create a basin object. A sediment object describes how sediment flux to the basin changes through time.
#' 
#' @param geometry an object of class [`geometry`].
#' @param startingVolume the starting sediment flux delivered at the left edge of the  
#'   model, in m * km.
#' @param netIncrease the amount by which sediment influx increases beginning of the  
#'   simulation to the end. Negative values indicate a decrease in sediment flux.
#' @param period the period of sediment flux cycles, in m.y.
#' @param amplitude the amplitude of sediment flux cycles, in m * km. The peak-to-peak  
#'   sediment flux change is twice the amplitude.
#' @param symmetry a value ranging from 0 to 1 describing the symmetry of the sediment  
#'   flux cycles. 0.5 is a classic cosine wave; smaller values have increasingly faster  
#'   initial change in rates and slower subsequent changes in rates; larger values  
#'   indicate the opposite.
#' @param phase a string with possible values of "falling", "rising", "highPoint", and  
#'   "lowPoint" describing the starting position of cyclic sediment flux.
#' @param shape a dimensionless value of 0 or greater that describes the shape of the  
#'   sediment flux cycle. 0 corresponds to a classic cosine curve, and larger values cause  
#'   the curve to become increasingly squared-off.
#' @param x,object an object of class `sediment`.
#' @param ... additional arguments to be passed.
#' 
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, 
#'   marginWidth=600, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' sedi <- sediment(geometry=geom, startingVolume=60)
#' print(sedi)
#' summary(sedi)
#' plot(sedi)
#' 
#' @rdname sediment
#' @export sediment
#' 
#' @return sediment returns an object of class "sediment", which includes print, summary, and plot methods.
#'
#' A sediment object is a list of two items. The first item (`parameters`) is a list of the arguments used in creating the sediment object. The second item (`timeSeries`) is a two-column data frame that records the time (`timePoint`, in m.y.) and the volume (`volume`, in m * km) of sediment supplied to the basin.
#' 


sediment <- function(geometry, startingVolume=60, netIncrease=0, period=1, amplitude=0, symmetry=0.5, phase=c('rising', 'falling', 'highPoint', 'lowPoint'), shape=0) {
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

#' @rdname sediment
#' @export

plot.sediment <- function(x, ...) {
	plot(x$timeSeries$timePoint, x$timeSeries$volume, type='l', las=1, xlab="model time (m.y.)", ylab="sediment flux", ...)
}

#' 
#' @rdname sediment
#' @export

print.sediment <- function(x, ...) {
	x
}

#' 
#' @rdname sediment
#' @export

summary.sediment <- function(object, ...) {
	cat("Starting volume of sediment (startingVolume):     ", object$parameters$startingVolume, "m^2*km\n")
	cat("Net increase in sediment (netIncrease):           ", object$parameters$netIncrease, "m^2*km\n")
	cat("Period of cyclical sediment input (period):       ", object$parameters$period, "m.y.\n")
	cat("Amplitude of cyclical sediment input (amplitude): ", object$parameters$amplitude, "m^2*km\n")
	cat("Symmetry of cyclical sediment input (symmetry):   ", object$parameters$symmetry, "(dimensionless, 0 to 1)\n")
	cat("Phase of cyclical sediment input (phase):         ", object$parameters$phase, "\n")
	cat("Shape of cyclical sediment input (shape):         ", object$parameters$shape, "(dimensionless, 0 to infinity)\n")
	cat("timePoints:                                       ", length(object$timeSeries$timePoint), "\n")
	cat("minimum sedment volume:                           ", min(object$timeSeries$volume), "m^2*km\n")
	cat("maximum sediment volume:                          ", max(object$timeSeries$volume), "m^2*km\n")
}
