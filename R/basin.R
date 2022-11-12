#' @title Basin Simulation
#'
#' @description Create an object of class "basin".
#'
#' @details Creates a basin object, which records the history of sedimentation, elevation,  
#'   and water depth across a sedimentary basin. Requires geometry, subsidence, eustasy,  
#'   and sediment objects. Sedimentary basins can take several minutes to run.
#' 
#' @param geometry an object of class [`geometry`].
#' @param subsidence an object of class [`subsidence`].
#' @param eustasy an object of class [`eustasy`].
#' @param sediment an object of class [`sediment`].
#' @param searchWidth a distance (in km) for the distance over which an optimal shore is  
#'   found. Generally should not be changed.
#' @param x,object an object of class `basin`.
#' @param setting a string ('valley' or 'interfluve') that specifies the location of a  
#'   basin cross-section.
#' @param addLegend a boolean indicating whether a legend for coastal plain and marine  
#'   facies should be added to the basin plot.
#' @param xlim,ylim Optional setting to set the range of the x-axis and y-axis that is  
#'   displayed.
#' @param ... additional arguments to be passed.

#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, 
#'   marginWidth=600, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' subs <- subsidence(geometry=geom, startingLeft=0.0, startingRight=1.0)
#' eust <- eustasy(geometry=geom, netRise=30.0)	
#' sedi <- sediment(geometry=geom, startingVolume=120, netIncrease=120)
#' sedBasin <- basin(geometry=geom, subsidence=subs, eustasy=eust, sediment=sedi)
#' summary(sedBasin)
#' print(sedBasin)
#' plot(sedBasin, setting="valley")
#' 
#' @rdname basin
#' @export basin
#' 
#' @return basin returns an object of class "basin", which includes print, summary, and plot methods.
#'
#' A basin object is a list consisting of 14 objects. The first element to the list, `parameters` contains a series of lists reflecting the inputs given by the geometry, subsidence, eustasy, and sediment objects. The remaining 13 items in the list are the outputs of the model. `positions` are the locations in the basin, in km from the left edge. `timePoints` are the model times that were simulated, in m.y. `elevationProfileValley` and `elevationProfileInterfluve` are matrices that record the elevation along those two profiles at every time step in the simulation; rows correspond to time points and columns correspond to positions in the basin. `sedimentAccumulatedValley` and `sedimentAccumulatedInterfluve` record the thickness of sediment accumulated (in meters) at every location and time point; rows correspond to time points and columns correspond to positions in the basin. `hiatusValley` and `hiatusInterfluve` are boolean matrices that record whether a hiatus was forming at any point at any time in the basin; rows correspond to time points and columns correspond to positions in the basin. `subsidenceRate` records the subsidence rate (in m / m.y.) at every location and time in the basin; rows correspond to time points and columns correspond to positions in the basin. `eustasy` is a vector of the position of sea level (in m) at each time point. `integratedSediment` is a vector of the total sediment volume deposited in the basin (in m * km) at each time point. `targetSediment` is a vector of the amount of sediment that was intended to be deposited at each time point; this should match `integratedSediment` closely. Finally, `shore` records the position of the shore at each time point, in km from the left edge of the basin.
#' 


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
		
		# Detecting the partitioning / near-zero sediment spikes
		if (bestFitValley$sedimentVolume / sediment$timeSeries$volume[i] < 0.5) {
			message(paste("insufficient sediment deposited at time point", i))
			message(paste("target volume:", sediment$timeSeries$volume[i], "; actual volume:", bestFitValley$sedimentVolume))
		}
	
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

#' @rdname basin
#' @export
plot.basin <- function(x, setting=c('valley', 'interfluve'), addLegend=TRUE, xlim=NULL, ylim=NULL, ...) {
	debugPlot=FALSE
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
	
	plot(x$positions, startSurface, type='n', xlim=xlim, ylim=ylim, xlab='Distance (km)', ylab='Elevation (m)', las=1, ...)

	# all sediment, filled with color
	sedimentsY <- c(startSurface, rev(endSurface))
	sedimentsX <- c(x$positions, rev(x$positions))
	graphics::polygon(sedimentsX, sedimentsY, col=marineColor)

	# draw a colored envelope around the coastal plain sediments
	startShoreIndex <- indexForPosition(positionKm=x$shore[1], marginWidth=x$parameters$marginWidth, deltaX=x$parameters$deltaX)
	endShoreIndex <- indexForPosition(positionKm=x$shore[endIndex], marginWidth=x$parameters$marginWidth, deltaX=x$parameters$deltaX)
	startCoastalPlain <- startSurface[1:startShoreIndex]
	endCoastalPlain <- endSurface[1:endShoreIndex]
	shoreCoords <- shores(basin=x, setting=setting)
	coastalPlainY <- c(startCoastalPlain, shoreCoords$y, rev(endCoastalPlain))
	coastalPlainX <- c(x$positions[1:startShoreIndex], shoreCoords$x, x$positions[endShoreIndex:1])
	graphics::polygon(coastalPlainX, coastalPlainY, col=coastalPlainColor, lwd=0.5)
	
	# add sea level at the end of the simulation
	endShoreX <- x$positions[endShoreIndex]
	graphics::segments(endShoreX, 0, x$parameters$marginWidth, 0, col=seaLevelColor)

	# optionally, add a legend for the colors
	if (addLegend) {
		plotTop <- max(ylim)
		plotWidth <- max(x$positions) - min(x$positions)
		boxTop <- 0.95 * plotTop
		boxLeft <- 0.75 * plotWidth
		graphics::legend(boxLeft, boxTop, fill=c(coastalPlainColor, marineColor), pt.cex=2, legend=c('coastal plain', 'marine'), box.col='white', bg='white')
	}
	
	if (debugPlot) {
		graphics::points(x$positions, startSurface, type='l', col='blue', lwd=2)
		graphics::points(x$positions, endSurface, type='l', col='red', lwd=2)
		shoreCoords <- shores(basin=x, setting='valley')
		graphics::points(shoreCoords$x, shoreCoords$y, type='l', col='yellow', lwd=2)
		shoreCoords <- shores(basin=x, setting='interfluve')
		graphics::points(shoreCoords$x, shoreCoords$y, type='l', col='yellow', lwd=2, lty='dotted')
	}
}

#' 
#' @rdname basin
#' @export
print.basin <- function(x, ...) {
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

#' 
#' @rdname basin
#' @export
summary.basin <- function(object, ...) {
	cat("Geometry\n")
	cat("fallLineX:         ", object$parameters$fallLineX, "km\n")
	cat("fallLineY:         ", object$parameters$fallLineY, "m\n")
	cat("shoreX:            ", object$parameters$shoreX, "km\n")
	cat("shoreY:            ", object$parameters$shoreY, "m\n")
	cat("deltaWidth:        ", object$parameters$deltaWidth, "km\n")
	cat("marginWidth:       ", object$parameters$marginWidth, "km\n")
	cat("deltaToeY:         ", object$parameters$deltaToeY, "m\n")
	cat("deltaX:            ", object$parameters$deltaX, "m\n")
	cat("nonMarAlpha:       ", object$parameters$nonMarAlpha, "(dimensionless)\n")
	cat("marineAlpha:       ", object$parameters$marineAlpha, "(dimensionless)\n")
	cat("duration:          ", object$parameters$duration, "m.y.\n")
	cat("timeStep:          ", object$parameters$timeStep, "m.y.\n")
	cat("\nSubsidence\n")
	cat("startingLeft:      ", object$parameters$subsidence$startingLeft, "m/m.y.\n")
	cat("startingRight:     ", object$parameters$subsidence$startingRight, "m/m.y.\n")
	cat("netChangeFactor:   ", object$parameters$subsidence$netChangeFactor, "(dimensionless)\n")
	cat("period:            ", object$parameters$subsidence$period, "m.y.\n")
	cat("amplitude:         ", object$parameters$subsidence$amplitude, "m\n")
	cat("symmetry:          ", object$parameters$subsidence$symmetry, "(dimensionless, 0 to 1\n")
	cat("phase:             ", object$parameters$subsidence$phase, "\n")
	cat("shape:             ", object$parameters$subsidence$shape, "(dimensionless, 0 to infinity)\n")
	cat("time steps:        ", nrow(object$subsidenceRate), "\n")
	cat("spatial positions: ", ncol(object$subsidenceRate), "\n")
	cat("minimum rate:      ", min(object$subsidenceRate), "m/m.y.\n")
	cat("maximum  rate:     ", max(object$subsidenceRate), "m/m.y.\n")
	cat("\nEustasy\n")
	cat("netRise:           ", object$parameters$eustasy$netRise, "m\n")
	cat("period):           ", object$parameters$eustasy$period, "m.y.\n")
	cat("amplitude:         ", object$parameters$eustasy$amplitude, "m\n")
	cat("symmetry:          ", object$parameters$eustasy$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:             ", object$parameters$eustasy$phase, "\n")
	cat("shape:             ", object$parameters$eustasy$shape, "(dimensionless, 0 to infinity)\n")
	cat("minimum sea level: ", min(object$eustasy), "m\n")
	cat("maximum sea level: ", max(object$eustasy), "m\n")
	cat("\nSediment\n")
	cat("startingVolume:    ", object$parameters$sedimentFluobject$startingVolume, "m^2*km\n")
	cat("netIncrease:       ", object$parameters$sedimentFluobject$netIncrease, "m^2*km\n")
	cat("period:            ", object$parameters$sedimentFluobject$period, "m.y.\n")
	cat("amplitude:         ", object$parameters$sedimentFluobject$amplitude, "m^2*km\n")
	cat("symmetry:          ", object$parameters$sedimentFluobject$symmetry, "(dimensionless, 0 to 1)\n")
	cat("phase:             ", object$parameters$sedimentFluobject$phase, "\n")
	cat("shape:             ", object$parameters$sedimentFluobject$shape, "(dimensionless, 0 to infinity)\n")
	cat("timePoints:        ", length(object$targetSediment), "\n")
	cat("minimum volume:    ", min(object$targetSediment), "m^2*km\n")
	cat("maximum  volume:   ", max(object$targetSediment), "m^2*km\n")
}


