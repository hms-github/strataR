#' Create a basin object
#'
#' Creates a basin object, which records the history of sedimentation, elevation, and water depth across a sedimentary basin. Requires geometry, subsidence, eustasy, and sediment objects. Sedimentary basins can take several minutes to run.
#' 
#' @title subsidence: The subsidence function
#' @param geometry a geometry object
#' @param subsidence a subsidence object
#' @param eustasy a eustasy object
#' @param sediment a sediment object
#' @param searchWidth a distance (in km) for the distance over which an optimal shore is found. Generally should not be changed.
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, marginWidth=500, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' subs <- subsidence(geometry=geom, startingLeft=0.0, startingRight=1.0)
#' eust <- eustasy(geometry=geom, netRise=30.0)	
#' sedi <- sediment(geometry=geom, startingVolume=120, netIncrease=120)
#' sedBasin <- basin(geometry=geom, subsidence=subs, eustasy=eust, sediment=sedi)
#' 
#' @rdname basin
#' @export basin

basin <- function(geometry, subsidence, eustasy, sediment, searchWidth=10) {
	timeStep <- geometry$timeStep
	eus <- eustasy$timeSeries
	shoreX <- geometry$shoreX
	nonMarAlpha <- geometry$nonMarAlpha
	marineAlpha <- geometry$marineAlpha
	
	integrityCheck(geometry, subsidence, eustasy, sediment)
		
	numTimePoints <- nrow(eus)                                # includes t0 already
	numSpatialPoints <- (geometry$marginWidth - geometry$fallLineX ) / geometry$deltaX + 1   # +1 to include x0
	positions <- seq(0, geometry$marginWidth, geometry$deltaX)

	elevationProfileValley <- matrix(nrow=numTimePoints, ncol=numSpatialPoints)
	elevationProfileInterfluve <- matrix(nrow=numTimePoints, ncol=numSpatialPoints)
	sedimentAccumulatedValley <- matrix(nrow=numTimePoints, ncol=numSpatialPoints)
	sedimentAccumulatedInterfluve <- matrix(nrow=numTimePoints, ncol=numSpatialPoints)
	hiatusValley <- matrix(data=FALSE, nrow=numTimePoints, ncol=numSpatialPoints)
	hiatusInterfluve <- matrix(data=FALSE, nrow=numTimePoints, ncol=numSpatialPoints)
	integratedSediment <- rep(NA, numTimePoints)
	shore <- rep(NA, numTimePoints)

	initial <- generateInitialProfile(fallLineX=geometry$fallLineX, fallLineY=geometry$fallLineY, shoreX=shoreX, shoreY=geometry$shoreY, deltaWidth=geometry$deltaWidth, marginWidth=geometry$marginWidth, deltaToeY=geometry$deltaToeY, nonMarAlpha=nonMarAlpha, marineAlpha=marineAlpha, deltaX=geometry$deltaX)

	elevationProfileValley[1, ] <- initial$y
	elevationProfileInterfluve[1, ] <- initial$y
	sedimentAccumulatedValley[1, ] <- rep(0, numSpatialPoints)
	sedimentAccumulatedInterfluve[1, ] <- rep(0, numSpatialPoints)
	integratedSediment[1] <- 0
	shore[1] <- shoreX

	for (i in 2:numTimePoints) {    # start at 2 because 1 is the initial profile
		shoreX <- shore[i-1]
		
		# deflect the previous profile
		deflectedValley <- elevationProfileValley[i-1, ]
		deflectedInterfluve <- elevationProfileInterfluve[i-1, ]
		eustaticRise <- eus$seaLevel[i] - eus$seaLevel[i-1]  # rise in sea level
		subside <- subsidence$rates[i, ] * timeStep
		deflectedValley <- deflectedValley + subside - eustaticRise
		deflectedInterfluve <- deflectedInterfluve + subside - eustaticRise
			
		fallLineY <- elevationProfileValley[i-1, 1] - eustaticRise

		# calculate the new profile, which will stay at grade
		deflected <- list(x=initial$x, y=deflectedValley)
		bestFitValley <- optimalProfile(deflected=deflected, previousShoreX=shoreX, deltaWidth=geometry$deltaWidth, marginWidth=geometry$marginWidth, deltaX=geometry$deltaX, searchWidth=searchWidth, targetSedVolume=sediment$timeSeries$volume[i], nonMarAlpha=nonMarAlpha, marineAlpha=marineAlpha, fallLineY=fallLineY)
		elevationProfileValley[i, ] <- bestFitValley$fittedProfile$y
		elevationProfileInterfluve[i, ] <- pmax(bestFitValley$fittedProfile$y, deflectedInterfluve)
	
		# calculate the sedimentation
		sedimentAccumulatedValley[i, ] <- elevationProfileValley[i, ] - deflectedValley
		sedimentAccumulatedInterfluve[i, ] <- elevationProfileInterfluve[i, ] - deflectedInterfluve
		
		# locate unconformities
		hiatusThreshold <- 0.01
		hiatusValley[i, ] <- (sedimentAccumulatedValley[i, ] <= hiatusThreshold)
		hiatusInterfluve[i, ] <- (sedimentAccumulatedInterfluve[i, ] <= hiatusThreshold)
		
		# save the total volume of sediment along the incised valley 
		integratedSediment[i] <- bestFitValley$sedimentVolume
		
		# save the current shore
		shore[i] <- bestFitValley$shore
	}
	
	positions <- seq(geometry$fallLineX, geometry$marginWidth, geometry$deltaX)
		
	parameters <- list(
		fallLineX=geometry$fallLineX, 
		fallLineY=geometry$fallLineY, 
		shoreX=geometry$shoreX, 
		shoreY=geometry$shoreY, 
		deltaWidth=geometry$deltaWidth, 
		marginWidth=geometry$marginWidth, 
		deltaToeY=geometry$deltaToeY, 
		deltaX=geometry$deltaX,
		duration=geometry$duration, 
		timeStep=geometry$timeStep, 
		nonMarAlpha=nonMarAlpha, 
		marineAlpha=marineAlpha, 
		subsidence=subsidence$parameters,
		eustasy=eustasy$parameters,
		sedimentFlux=sediment$parameters
	)
	
	results <- list(
		parameters=parameters, 
		positions=positions, 
		timePoints=eus$timePoint, 
		elevationProfileValley=elevationProfileValley,  
		elevationProfileInterfluve=elevationProfileInterfluve, 
		sedimentAccumulatedValley=sedimentAccumulatedValley,  
		sedimentAccumulatedInterfluve=sedimentAccumulatedInterfluve, 
		hiatusValley=hiatusValley, 
		hiatusInterfluve=hiatusInterfluve, 
		subsidenceRate=subsidence$rates,
		eustasy=eus$seaLevel, 
		integratedSediment=integratedSediment, 
		targetSediment=sediment$timeSeries$volume,
		shore=shore
	)
			
	class(results) <- "basin"
	return(results)
}

#' @return \code{NULL}
#' 
#' @rdname basin
#' @export
plot.basin <- function(x, setting=c('valley', 'interfluve'), addLegend=TRUE, xlim=NULL, ylim=NULL, debugPlot=FALSE) {
	setting <- match.arg(setting)
	endIndex <- length(x$timePoints)
	endTime <- x$timePoints[endIndex]
	
	marineColor <- 'tan'
	coastalPlainColor <- 'olivedrab3'
	seaLevelColor <- 'dodgerblue'
	
	deflection <- finalDeflection(x)
	
	elevation <- 1
	if (setting == 'valley') {
	 	elevation <- x$elevationProfileValley
	} else if (setting == 'interfluve') {
		elevation <- x$elevationProfileInterfluve
	} else {
	     warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
	}
	
	startSurface <- elevation[1, ] - deflection[1, ]
	endSurface <- elevation[endIndex, ]
	
	# handle any erosion of the starting surface at a valley
	if (setting == 'valley') {
		for (timeIndex in 2:endIndex) {
			overlyingSurface <- elevation[timeIndex, ]
			deflectedOverlyingSurface <- overlyingSurface - deflection[timeIndex, ] 
			startSurface <- pmin(startSurface, deflectedOverlyingSurface)
		}
	}

	# If no specific range is set for x and y, create the default plot limits
	if (is.null(xlim)) {
		xlim=range(x$positions)
	}
	if (is.null(ylim)) {
		ylim <- c(min(endSurface, startSurface), max(endSurface, startSurface))
	}
	
	plot(x$positions, startSurface, type='n', xlim=xlim, ylim=ylim, xlab='Distance (km)', ylab='Elevation (m)', las=1)

	# all sediment, filled with color
	sedimentsY <- c(startSurface, rev(endSurface))
	sedimentsX <- c(x$positions, rev(x$positions))
	polygon(sedimentsX, sedimentsY, col=marineColor)

	# draw a colored envelope around the coastal plain sediments
	startShoreIndex <- indexForPosition(positionKm=x$shore[1], marginWidth=x$parameters$marginWidth, deltaX=x$parameters$deltaX)
	endShoreIndex <- indexForPosition(positionKm=x$shore[endIndex], marginWidth=x$parameters$marginWidth, deltaX=x$parameters$deltaX)
	startCoastalPlain <- startSurface[1:startShoreIndex]
	endCoastalPlain <- endSurface[1:endShoreIndex]
	shoreCoords <- shores(basin=x, setting=setting)
	coastalPlainY <- c(startCoastalPlain, shoreCoords$y, rev(endCoastalPlain))
	coastalPlainX <- c(x$positions[1:startShoreIndex], shoreCoords$x, x$positions[endShoreIndex:1])
	polygon(coastalPlainX, coastalPlainY, col=coastalPlainColor, lwd=0.5)
	
	# add sea level at the end of the simulation
	endShoreX <- x$positions[endShoreIndex]
	segments(endShoreX, 0, x$parameters$marginWidth, 0, col=seaLevelColor)

	# optionally, add a legend for the colors
	if (addLegend) {
		plotTop <- max(ylim)
		plotWidth <- max(x$positions) - min(x$positions)
		boxTop <- 0.95 * plotTop
		boxLeft <- 0.75 * plotWidth
		legend(boxLeft, boxTop, fill=c(coastalPlainColor, marineColor), pt.cex=2, legend=c('coastal plain', 'marine'), box.col='white', bg='white')
	}
	
	if (debugPlot) {
		points(x$positions, startSurface, type='l', col='blue', lwd=2)
		points(x$positions, endSurface, type='l', col='red', lwd=2)
		shoreCoords <- shores(basin=x, setting='valley')
		points(shoreCoords$x, shoreCoords$y, type='l', col='yellow', lwd=2)
		shoreCoords <- shores(basin=x, setting='interfluve')
		points(shoreCoords$x, shoreCoords$y, type='l', col='yellow', lwd=2, lty='dotted')
	}
}

#' @return \code{NULL}
#' 
#' @rdname basin
#' @export
print.basin <- function(x) {
	cat("Geometry\n")
	cat("fallLineX:         ", x$parameters$fallLineX, "km\n")
	cat("fallLineY:         ", x$parameters$fallLineY, "m\n")
	cat("shoreX:            ", x$parameters$shoreX, "km\n")
	cat("shoreY:            ", x$parameters$shoreY, "m\n")
	cat("deltaWidth:        ", x$parameters$deltaWidth, "km\n")
	cat("marginWidth:       ", x$parameters$marginWidth, "km\n")
	cat("deltaToeY:         ", x$parameters$deltaToeY, "m\n")
	cat("deltaX:            ", x$parameters$deltaX, "m\n")
	cat("nonMarAlpha:       ", x$parameters$nonMarAlpha, "(dimensionless)\n")
	cat("marineAlpha:       ", x$parameters$marineAlpha, "(dimensionless)\n")
	cat("duration:          ", x$parameters$duration, "m.y.\n")
	cat("timeStep:          ", x$parameters$timeStep, "m.y.\n")
	cat("\nSubsidence\n")
	cat("startingLeft:      ", x$parameters$subsidence$startingLeft, "m/m.y.\n")
	cat("startingRight:     ", x$parameters$subsidence$startingRight, "m/m.y.\n")
	cat("netChangeFactor:   ", x$parameters$subsidence$netChangeFactor, "(dimensionless)\n")
	cat("period:            ", x$parameters$subsidence$period, "m.y.\n")
	cat("amplitude:         ", x$parameters$subsidence$amplitude, "m\n")
	cat("symmetry:          ", x$parameters$subsidence$symmetry, "(dimensionless, 0 to 1\n")
	cat("phase:             ", x$parameters$subsidence$phase, "\n")
	cat("shape:             ", x$parameters$subsidence$shape, "(dimensionless, 0 to infinity)\n")
	cat("time steps:        ", nrow(x$subsidenceRate), "\n")
	cat("spatial positions: ", ncol(x$subsidenceRate), "\n")
	cat("minimum rate:      ", min(x$subsidenceRate), "m/m.y.\n")
	cat("maximum  rate:     ", max(x$subsidenceRate), "m/m.y.\n")
	cat("\nEustasy\n")
	cat("netRise:           ", x$parameters$eustasy$netRise, "m\n")
	cat("period):           ", x$parameters$eustasy$period, "m.y.\n")
	cat("amplitude:         ", x$parameters$eustasy$amplitude, "m\n")
	cat("symmetry:          ", x$parameters$eustasy$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:             ", x$parameters$eustasy$phase, "\n")
	cat("shape:             ", x$parameters$eustasy$shape, "(dimensionless, 0 to infinity)\n")
	cat("minimum sea level: ", min(x$eustasy), "m\n")
	cat("maximum sea level: ", max(x$eustasy), "m\n")
	cat("\nSediment\n")
	cat("startingVolume:    ", x$parameters$sedimentFlux$startingVolume, "m^2*km\n")
	cat("netIncrease:       ", x$parameters$sedimentFlux$netIncrease, "m^2*km\n")
	cat("period:            ", x$parameters$sedimentFlux$period, "m.y.\n")
	cat("amplitude:         ", x$parameters$sedimentFlux$amplitude, "m^2*km\n")
	cat("symmetry:          ", x$parameters$sedimentFlux$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:             ", x$parameters$sedimentFlux$phase, "\n")
	cat("shape:             ", x$parameters$sedimentFlux$shape, "(dimensionless, 0 to infinity)\n")
	cat("timePoints:        ", length(x$targetSediment), "\n")
	cat("minimum volume:    ", min(x$targetSediment), "m^2*km\n")
	cat("maximum  volume:   ", max(x$targetSediment), "m^2*km\n")
}

#' @return \code{NULL}
#' 
#' @rdname basin
#' @export
summary.basin <- function(x) {
	cat("Geometry\n")
	cat("fallLineX:         ", x$parameters$fallLineX, "km\n")
	cat("fallLineY:         ", x$parameters$fallLineY, "m\n")
	cat("shoreX:            ", x$parameters$shoreX, "km\n")
	cat("shoreY:            ", x$parameters$shoreY, "m\n")
	cat("deltaWidth:        ", x$parameters$deltaWidth, "km\n")
	cat("marginWidth:       ", x$parameters$marginWidth, "km\n")
	cat("deltaToeY:         ", x$parameters$deltaToeY, "m\n")
	cat("deltaX:            ", x$parameters$deltaX, "m\n")
	cat("nonMarAlpha:       ", x$parameters$nonMarAlpha, "(dimensionless)\n")
	cat("marineAlpha:       ", x$parameters$marineAlpha, "(dimensionless)\n")
	cat("duration:          ", x$parameters$duration, "m.y.\n")
	cat("timeStep:          ", x$parameters$timeStep, "m.y.\n")
	cat("\nSubsidence\n")
	cat("startingLeft:      ", x$parameters$subsidence$startingLeft, "m/m.y.\n")
	cat("startingRight:     ", x$parameters$subsidence$startingRight, "m/m.y.\n")
	cat("netChangeFactor:   ", x$parameters$subsidence$netChangeFactor, "(dimensionless)\n")
	cat("period:            ", x$parameters$subsidence$period, "m.y.\n")
	cat("amplitude:         ", x$parameters$subsidence$amplitude, "m\n")
	cat("symmetry:          ", x$parameters$subsidence$symmetry, "(dimensionless, 0 to 1\n")
	cat("phase:             ", x$parameters$subsidence$phase, "\n")
	cat("shape:             ", x$parameters$subsidence$shape, "(dimensionless, 0 to infinity)\n")
	cat("time steps:        ", nrow(x$subsidenceRate), "\n")
	cat("spatial positions: ", ncol(x$subsidenceRate), "\n")
	cat("minimum rate:      ", min(x$subsidenceRate), "m/m.y.\n")
	cat("maximum  rate:     ", max(x$subsidenceRate), "m/m.y.\n")
	cat("\nEustasy\n")
	cat("netRise:           ", x$parameters$eustasy$netRise, "m\n")
	cat("period):           ", x$parameters$eustasy$period, "m.y.\n")
	cat("amplitude:         ", x$parameters$eustasy$amplitude, "m\n")
	cat("symmetry:          ", x$parameters$eustasy$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:             ", x$parameters$eustasy$phase, "\n")
	cat("shape:             ", x$parameters$eustasy$shape, "(dimensionless, 0 to infinity)\n")
	cat("minimum sea level: ", min(x$eustasy), "m\n")
	cat("maximum sea level: ", max(x$eustasy), "m\n")
	cat("\nSediment\n")
	cat("startingVolume:    ", x$parameters$sedimentFlux$startingVolume, "m^2*km\n")
	cat("netIncrease:       ", x$parameters$sedimentFlux$netIncrease, "m^2*km\n")
	cat("period:            ", x$parameters$sedimentFlux$period, "m.y.\n")
	cat("amplitude:         ", x$parameters$sedimentFlux$amplitude, "m^2*km\n")
	cat("symmetry:          ", x$parameters$sedimentFlux$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:             ", x$parameters$sedimentFlux$phase, "\n")
	cat("shape:             ", x$parameters$sedimentFlux$shape, "(dimensionless, 0 to infinity)\n")
	cat("timePoints:        ", length(x$targetSediment), "\n")
	cat("minimum volume:    ", min(x$targetSediment), "m^2*km\n")
	cat("maximum  volume:   ", max(x$targetSediment), "m^2*km\n")
}


