#' Create a subsidence object
#'
#' Creates a subsidence object, which is needed to create a basin object. A subsidence object describes how subsidence changes across the sedimentary basin and through time
#' 
#' @title subsidence: The subsidence function
#' @param geometry a geometry object
#' @param startingLeft the starting subsidence rate at the left edge of the model, where sediment is introduced, in m/m.y.
#' @param startingRight the starting subsidence rate at the right edge of the model, farthest from where sediment is introduced, in m/m.y.
#' @param netChangeFactor the factor by which subsidence changes from the beginning of the simulation to the end. For example, 0.5 would indicate a halving of subsidence rates, 1 would indicate no change in rates, and 2 would indicate a doubling.
#' @param period the period of subsidence cycles, in m.y.
#' @param amplitude the amplitude of subsidence cycles, in m. The peak-to-peak subsidence change is twice the amplitude.
#' @param symmetry a value ranging from 0 to 1 describing the symmetry of the subsidence cycles. 0.5 is a classic cosine wave; smaller values have increasingly faster initial change in rates and slower subsequent changes in rates; larger values indicate the opposite.
#' @param phase a string with possible values of rising, falling, highPoint, and lowPoint describing the starting position of cyclic subsidence.
#' @param shape a value of 0 or greater that describes the shape of the subsidence cycle. 0 corresponds to a classic cosine curve, and larger values cause the curve to become increasingly squared-off.
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, marginWidth=500, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' subs <- subsidence(geometry=geom, startingLeft=0.0, startingRight=1.0)
#' 
#' @rdname subsidence
#' @export subsidence

subsidence <- function(geometry, startingLeft=0.0, startingRight, netChangeFactor=1, period=1, amplitude=0, symmetry=0.5, phase=c('rising', 'falling', 'highPoint', 'lowPoint'), shape=0) {
	# geometry is the object supplied by setGeometry()
	# all times (period, duration, timeStep, timePoint) are in m.y.
	# all rates (startingLeft, startingRight) are in m/m.y. and reflect the subsidence rates at the left edge (at the sediment source) and right edge of the simulation
	# netIncreaseFactor is multiplier expressing how much the final subsidence rate (anywhere) is greater than the initial subsidence rate at that spot, as per m.y.; 0 means no increase, 1 is 100% increase, etc.
	# phase is one of four values - "falling", "rising", "highPoint", "lowPoint" - that describe whether the sine wave starts on the falling inflection point, rising inflection point, highest point, or lowest point on the sine wave.
	# shape is dimensionless, 0 is a sine wave, higher values are flattened sine waves
	
	if (netChangeFactor <= 0) {
		stop("Argument netChangeFactor must have a value greater than zero", call.=FALSE)
	}
	
	phase <- match.arg(phase)
	duration <- geometry$duration
	timeStep <- geometry$timeStep
	timePoint <- seq(0, duration, timeStep)
	initialSubsRates <- crossShelfSubsidence(subsRateLeft=startingLeft, subsRateRight=startingRight, geometry=geometry)$rates
	
	sineComponentLeft <- 0
	sineComponentRight <- 0
	if (startingLeft > startingRight) {  # foreland
		sineComponentLeft <- flexSin(timePoint, period=period, amplitude=amplitude, symmetry=symmetry, phase=phase, shape=shape)
	} else {                             # passive margin
		sineComponentRight <- flexSin(timePoint, period=period, amplitude=amplitude, symmetry=symmetry, phase=phase, shape=shape)
	}
	
	netChangeFactor <- netChangeFactor - 1 # so that an input of 1 causes no change
		 
	trendLeft <- startingLeft + startingLeft * netChangeFactor * timePoint/duration
	subsLeft <- sineComponentLeft + trendLeft
	
	trendRight <- startingRight + startingRight * netChangeFactor * timePoint/duration
	subsRight <- sineComponentRight + trendRight
	
	subsidence <- matrix(0, nrow=length(timePoint), ncol=length(initialSubsRates))
	subsidence[1, ] <- initialSubsRates
	for (i in 1:length(timePoint)) {
		subsTimePoint <- crossShelfSubsidence(subsRateLeft=subsLeft[i], subsRateRight=subsRight[i], geometry)
		subsidence[i, ] <- subsTimePoint$rates
	}
	
	parameters <- list(startingLeft=startingLeft, startingRight=startingRight, netIncreaseFactor=netChangeFactor, period=period, amplitude=amplitude, symmetry=symmetry, phase=phase, shape=shape)
	results <- list(parameters=parameters, rates=subsidence)
	class(results) <- "subsidence"
	return(results)
}

#' @return \code{NULL}
#' 
#' @rdname subsidence
#' @export

plot.subsidence <- function(x, geometry, type=c('lines', 'filled')) {
	type <- match.arg(type)
	numPositions <- ncol(x$rates)
	numTimeSteps <- nrow(x$rates)
	xaxis <- (0:(numPositions-1))*geometry$deltaX
	yaxis <- (0:(numTimeSteps-1))*geometry$timeStep
	if (type=='lines') {
		contour(xaxis, yaxis, -t(x$rates), main='subsidence rates (m / m.y.)', xlab='position (km)', ylab='model time (m.y.)', las=1)
	} else if (type=='filled') {
		filled.contour(xaxis, yaxis, -t(x$rates), main='subsidence rates (m / m.y.)', xlab='position (km)', ylab='model time (m.y.)', las=1)
	} else {
		warning("Invalid type, must be 'lines' or 'filled'", call.=FALSE, immediate.=TRUE)
	}
}

#' @return \code{NULL}
#' 
#' @rdname subsidence
#' @export

print.subsidence <- function(x) {
	cat("startingLeft:                ", x$parameters$startingLeft, "m/m.y.\n")
	cat("startingRight:               ", x$parameters$startingRight, "m/m.y.\n")
	cat("netChangeFactor:             ", x$parameters$netChangeFactor, "(dimensionless)\n")
	cat("period:                      ", x$parameters$period, "m.y.\n")
	cat("amplitude:                   ", x$parameters$amplitude, "m\n")
	cat("symmetry:                    ", x$parameters$symmetry, "(dimensionless, 0 to 1\n")
	cat("phase:                       ", x$parameters$phase, "\n")
	cat("shape:                       ", x$parameters$shape, "(dimensionless, 0 to infinity)\n")
	cat("Number of time steps:        ", nrow(x$rates), "\n")
	cat("Number of spatial positions: ", ncol(x$rates), "\n")
	cat("Minimum subsidence rate:     ", min(x$rates), "m/m.y.\n")
	cat("Maximum subsidence rate:     ", max(x$rates), "m/m.y.\n")
}

#' @return \code{NULL}
#' 
#' @rdname subsidence
#' @export

summary.subsidence <- function(x) {
	cat("Initial subsidence rate at left edge (startingLeft):                  ", x$parameters$startingLeft, "m/m.y.\n")
	cat("Initial subsidence rate at right edge (startingRight):                ", x$parameters$startingRight, "m/m.y.\n")
	cat("Factor of net change in subsidence rates over time (netChangeFactor): ", x$parameters$netChangeFactor, "(dimensionless)\n")
	cat("Period of cyclical subsidence (period):                               ", x$parameters$period, "m.y.\n")
	cat("Amplitude of cyclical subsidence (amplitude):                         ", x$parameters$amplitude, "m\n")
	cat("Symmetry of cyclical subsidence (symmetry):                           ", x$parameters$symmetry, "(dimensionless, 0 to 1)\n")
	cat("Phase of cyclical subsidence (phase):                                 ", x$parameters$phase, "\n")
	cat("Shape of cyclical subsidence (shape):                                 ", x$parameters$shape, "(dimensionless, 0 to infinity)\n")
}
