# ---------------------------------------------------------------------------------------
## UTILITY FUNCTIONS

# Obtain the x-index for a position within the basin
indexForPosition <- function(positionKm, marginWidth, deltaX) {
	# positionKm, marginWidth, and deltaX all in km
	allPositions <- seq(0, marginWidth, deltaX)
	deviations <- abs(allPositions - positionKm)
	index <- which(deviations == min(deviations))
	index <- index[1]   # in case of a tie, choose the left-most position index
	return(index)
}

# Calculate the amount each elevation profile has been deflected by subsidence and eustasy
# Returns a matrix of deflections with rows being time points (1 is oldest) and columns being positions (1 is leftmost in basin)
finalDeflection <- function(basin) {
	# reverse order of rows in subsidence, because we want to sum subsidence going backwards in time
	subsRate <- basin$subsidenceRate[nrow(basin$subsidenceRate):1, ]
	subsidence <- apply(subsRate, MARGIN=2, FUN=cumsum) * basin$parameters$timeStep
	# reverse it again so that row 1 is again the first time step,
	# multiply by -1 so that positive is a downward deflection
	subsidence <- -subsidence[nrow(subsidence):1, ]
	netEustaticRise <- basin$eustasy[length(basin$eustasy)] - basin$eustasy
	deflection <- subsidence + netEustaticRise
	deflection
}

# Return a vector of all elevations through time at a location in the basin
# Used by generateColumnFromBasin() in sedimentationModule
# Used by elevationTimePlot()
elevationAtLocation <- function(basin, locationKm, setting=c('valley', 'interfluve')) {
	setting <- match.arg(setting)
	locations <- seq(0, basin$parameter$marginWidth, basin$parameters$deltaX)

	# find the closest point, and in the case of a tie, get the first match
	locationIndex <- which(abs(locations-locationKm) == min(abs(locations-locationKm)))[1]

	allElevations <- 1
	if (setting == 'valley') {
	 	allElevations <- basin$elevationProfileValley
	} else if (setting == 'interfluve') {
		allElevations <- basin$elevationProfileInterfluve
	} else {
	     warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
	}

	elevation <- allElevations[ , locationIndex]
	return(elevation)
}

# Return a vector of sediment accumulation and erosion through time at a location in the basin
# Used by generateColumnFromBasin() in sedimentationModule
sedimentationAtLocation <- function(basin, locationKm, setting=c('valley', 'interfluve')) {
	setting <- match.arg(setting)
	locations <- seq(0, basin$parameter$marginWidth, basin$parameters$deltaX)

	# find the closest point, and in the case of a tie, get the first match
	locationIndex <- which(abs(locations-locationKm) == min(abs(locations-locationKm)))[1]

	allSedimentation <- 1
	if (setting == 'valley') {
	 	allSedimentation <- basin$sedimentAccumulatedValley
	} else if (setting == 'interfluve') {
		allSedimentation <- basin$sedimentAccumulatedInterfluve
	} else {
	     warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
	}

	sedimentation <- allSedimentation[ , locationIndex]
	return(sedimentation)
}

# Calculate total amounts of sediment in basin, amount in nonmarine settings, amount in marine settings, and the proportion of sediment that is nonmarine
# Used by partitioninPlot()
sedimentPartitioning <- function(basin) {
	# always calculated for valley, since that is where sedimentation is continuous
	allProfiles <- basin$elevationProfileValley
	timeStep <- basin$parameters$timeStep
	leftEdge <- 1
	rightEdge <- ncol(allProfiles)
	numTimes <- nrow(allProfiles)
	nonmarineVolume <- rep(0, numTimes)
	marineVolume <- rep(0, numTimes)
	totalVolume <- rep(0, numTimes)
	fractionNonmarine <- rep(0, numTimes)

	for (timeIndex in 2:numTimes) {
		lowerProfile <- allProfiles[timeIndex-1, ]
		lowerShoreIndex <- indexForPosition(positionKm=basin$shore[timeIndex-1], marginWidth=basin$parameters$marginWidth, deltaX=basin$parameters$deltaX)

		eustaticRise <- basin$eustasy[timeIndex] - basin$eustasy[timeIndex-1]  # rise in sea level
		deflectedLower <- lowerProfile
		subside <- basin$subsidenceRate[timeIndex, ] * timeStep
		deflectedLower <- deflectedLower + subside - eustaticRise

		upperProfile <- allProfiles[timeIndex, ]
		upperShoreIndex <- indexForPosition(positionKm=basin$shore[timeIndex], marginWidth=basin$parameters$marginWidth, deltaX=basin$parameters$deltaX)

		sediment <- upperProfile - deflectedLower
		sediment[sediment < 0] <- 0   # correct for any truncation

		if (lowerShoreIndex == upperShoreIndex) {
			nonmarineVolume[timeIndex] <- sum(sediment[leftEdge:upperShoreIndex])
			marineVolume[timeIndex] <- sum(sediment[(upperShoreIndex+1):rightEdge])
			totalVolume[timeIndex] <- sum(sediment)
		} else {
			alwaysNonmarine <- leftEdge : min(lowerShoreIndex, upperShoreIndex)
			alwaysMarine <- (max(lowerShoreIndex, upperShoreIndex) + 1) : rightEdge
			transitional <- (max(alwaysNonmarine)+1):(min(alwaysMarine-1))
			nonmarineVolume[timeIndex] <- sum(sediment[alwaysNonmarine]) + sum(sediment[transitional])/2
			marineVolume[timeIndex] <- sum(sediment[alwaysMarine]) + sum(sediment[transitional])/2
			totalVolume[timeIndex] <- sum(sediment)
		}
	}

	nonmarineVolume <- nonmarineVolume * basin$parameters$deltaX
	marineVolume <- marineVolume * basin$parameters$deltaX
	totalVolume <- totalVolume * basin$parameters$deltaX
	fractionNonmarine <- nonmarineVolume / totalVolume

	results <- data.frame(nonmarineVolume, marineVolume, totalVolume, fractionNonmarine)
	results <- results[-1, ]   # remove t0
	return(results)
}

# Return a vector of accommodation rates at the shoreline (in the valley) through time,
#   and optionally, a plot
shoreAccommodationRates <- function(basin) {
	elevation <- basin$elevationProfileValley
	timePoints <- basin$timePoints
	seaLevel <- basin$eustasy
	shoreX <- basin$shore

	numTimePoints <- length(timePoints)
	shoreMovement <- rep(0, numTimePoints)
	eustaticRise <- rep(0, numTimePoints)
	shoreSubsidence <- rep(0, numTimePoints)
	accommodation <- rep(0, numTimePoints)

	# start at t1 (timeIndex2), not t0 (no subsidence or eustatic change at t0)
	for (timeIndex in 2:numTimePoints) {
		previousShoreIndex <- indexForPosition(positionKm=shoreX[timeIndex-1], marginWidth=basin$parameters$marginWidth, deltaX=basin$parameters$deltaX)
		currentShoreIndex <- indexForPosition(positionKm=shoreX[timeIndex], marginWidth=basin$parameters$marginWidth, deltaX=basin$parameters$deltaX)
		shoreMovement[timeIndex] <- shoreX[timeIndex] - shoreX[timeIndex-1]

		# find the rate of subsidence at each of those shorelines
		# flip sign so positive shoreSubsidence means subsidence, negative is uplift
		shoreSubsidence[timeIndex] <- -(basin$subsidenceRate[timeIndex, currentShoreIndex] * basin$parameters$timeStep)

		# find the eustatic rise at each time step
		eustaticRise[timeIndex] <- seaLevel[timeIndex] - seaLevel[timeIndex-1]

		# calculate the rate of accommodation at the shore at each time step
		accommodation[timeIndex] <- shoreSubsidence[timeIndex] + eustaticRise[timeIndex]
	}

	# remove t0
	results <- data.frame(modelTime=timePoints[-1], shoreKm=shoreX[-1], shoreMovement=shoreMovement[-1], shoreSubsidence=shoreSubsidence[-1], eustaticRise=eustaticRise[-1], accommodation=accommodation[-1])

	return(results)
}

# Plots the accommodation rates at the shore, based on output of shoreAccommodationRates()
shoreAccommodationPlot <- function(basin, shoreRates) {
	shoreRates <- shoreAccommodationRates(basin)
	minRate <- min(0, shoreRates$accommodation)
	maxRate <- max(0, shoreRates$accommodation)
	plot(shoreRates$modelTime, shoreRates$accommodation, type='l', xlab='Model time (m.y.)', ylab='Accommodation rate (m/m.y.)', ylim=c(minRate, maxRate), las=1, lwd=2)

	# Zero-rate line
	graphics::abline(h=0, col='black')

	# Subsidence rates at the shore
	graphics::points(shoreRates$modelTime, shoreRates$shoreSubsidence, type='l', col='gray', lty='dotted')

	# Times in which the shore reverses its direction of movement
	# Note that this only finds times when the shore position is not changing
	# If it changes immediately from regression to transgression, or vice versa, without pausing for a time of no change, this will not find that reversal
	shoreMoveSign <- shoreRates$shoreMovement / abs(shoreRates$shoreMovement)
	shoreMoveSign[is.nan(shoreMoveSign)] <- 0
	graphics::abline(v=which(shoreMoveSign == 0)*basin$parameters$timeStep, col='gray', lty='dashed')
}

# Obtain the final x (in km) and y (elevation, in m) coordinates of all shores
# Used in basinPlot()
shores <- function(basin, setting=c('valley', 'interfluve')) {
	setting <- match.arg(setting)

	elevation <- 1
	deflection <- finalDeflection(basin)

	if (setting == 'valley') {
	 	elevation <- basin$elevationProfileValley
	} else if (setting == 'interfluve') {
		elevation <- basin$elevationProfileInterfluve
	} else {
	     warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
	}

	endIndex <- length(basin$timePoints)
	endTime <- basin$timePoints[endIndex]

	indices <- 1:length(basin$timePoints)
	shoreX <- basin$shore
	shoreY <- rep(0, length(indices))

	# find all the shorelines, lowest profile first
	timeIndex <- 1
	currentTime <- basin$timePoints[timeIndex]
	currentProfile <- elevation[timeIndex, ]
	currentShoreIndex <- indexForPosition(positionKm=shoreX[timeIndex], marginWidth=basin$parameters$marginWidth, deltaX=basin$parameters$deltaX)
	deflectedProfile <- currentProfile - deflection[timeIndex, ]
	shoreY[timeIndex] <- deflectedProfile[currentShoreIndex]

	for (timeIndex in 2:max(indices)) {
		currentTime <- basin$timePoints[timeIndex]
		currentProfile <- elevation[timeIndex, ]
		currentShoreIndex <- indexForPosition(positionKm=shoreX[timeIndex], marginWidth=basin$parameters$marginWidth, deltaX=basin$parameters$deltaX)
		deflectedProfile <- currentProfile - deflection[timeIndex, ]
		shoreY[timeIndex] <- deflectedProfile[currentShoreIndex]
		# in the valley where erosion can occur, correct all older shore positions to be no higher than the current profile
		if (setting == 'valley') {
			for (olderTimeIndex in 1:(timeIndex-1)) {
				previousShoreXindex <- indexForPosition(positionKm=shoreX[olderTimeIndex], marginWidth=basin$parameters$marginWidth, deltaX=basin$parameters$deltaX)
				if (shoreY[olderTimeIndex] >= deflectedProfile[previousShoreXindex]) {
					 shoreY[olderTimeIndex] <- deflectedProfile[previousShoreXindex]
				}
			}
		}

	}

	shoreCoordinates <- data.frame(x=shoreX, y=shoreY)
	shoreCoordinates
}




# ---------------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# Functions not generally called by user but by other functions

# Called by flexSin() to allow for asymmetry in the sine wave. Shouldn't be necessary for a user to call this.
symmetryScaling <- function(x, symmetry=0.5) {
	midpoint <- 0.5
	if (x <= symmetry) {
		x <- x / symmetry * midpoint
	} else {
		x <- (x - symmetry)/(1.0 - symmetry) * (1 - midpoint) + midpoint
	}
	return(x)
}

# Called by flexSin() to allow for a phase shift in the sine wave. Shouldn't be necessary for a user to call this.
phaseScaling <- function(x, phase=c("falling", "rising", "highPoint", "lowPoint"), symmetry=0.5) {
	if (phase=="highPoint") {       # a valid string was supplied
		phaseShift <- 0
	} else if (phase=="falling") {
		phaseShift <- symmetry - symmetry/2
	} else if (phase=="lowPoint") {
		phaseShift <- symmetry
	} else if (phase=="rising") {
		phaseShift <- symmetry + (1 - symmetry)/2
	} else  {                        # not a valid string
		stop("Invalid phase specified: must be \"falling\", \"rising\", \"highPoint\", or \"lowPoint\".", call.=FALSE)
	}
	xPhased <- (x + phaseShift)%%1
	return(xPhased)
}

# Used by eustasy(), subsidence(), and sediment() to create sinusoidal like waves that are squared off. This has now been subsumed into flexSin() which also handles asymmetry.
flatSin <- function(x, period=1, amplitude=1, phase=c("falling", "rising", "highPoint", "lowPoint"), shape=0) {
	phase <- phase[1]

	if (phase=="falling") {                      # a valid string was supplied
		phase <- period/2
	} else if (phase=="rising") {
		phase <- 0
	} else if (phase=="highPoint") {
		phase <- period/4
	} else if (phase=="lowPoint") {
		phase <- 3*period/4
	} else {
		if (is.numeric(phase) == TRUE) {         # a numeric value was supplied
			phase <- phase
		} else {                                 # neither a numeric value or valid string
			stop("Invalid phase specified: must be a number or \"falling\", \"rising\", \"highPoint\", or \"lowPoint\".", call.=FALSE)
		}
	}

	y <- amplitude * sqrt((1+shape^2) / (1+shape^2*sin((x+phase)/period*2*pi)^2)) * sin((x+phase)/period*2*pi)

	# shift curve so that it starts at zero
	y <- y - y[1]
	y
}

# Used as an input to fillBasin() and used by finalDeflection(), itself used in basinPlot()
crossShelfSubsidence <- function(subsRateLeft, subsRateRight, geometry) {
	# subsRateLeft and subsRateRight in m / m.y.
	positions <- seq(0, geometry$marginWidth, geometry$deltaX)
	rates <- (subsRateLeft-subsRateRight) / geometry$marginWidth * positions - subsRateLeft
	parameters <- list(subsRateLeft=subsRateLeft, subsRateRight=subsRateRight)
	results <- list(parameters=parameters, rates=rates)
	return(results)
}

# Produces a Caputo curve
# Curvature determine by alpha: 0.5 is slightly curved (~linear), 5.0 is strongly curved
# Paola uses 0.5
# Modified from what is in fluvialProfileModule.R to remove the needless shape argument
# Used by generateInitialProfile() and profileForShore(), ultimately by fillBasin()
height <- function(x, alpha) {
	shape <- 1  # required by dgamma, but seems to have no effect
	2 / stats::dgamma(alpha + 2, shape) * (1 - x) ^ (alpha + 1)
}

# Rescale a Caputo curve in the x-direction
# Values monotonically increase, as in with distance from fluvial source
# Used by generateInitialProfile() and profileForShore(), ultimately by fillBasin()
xRescale <- function(x, leftX, rightX) {
	maxX <- max(x)
	minX <- min(x)
	rescaled <- x * (rightX - leftX) / (maxX - minX) + leftX
	rescaled
}

# Rescale a Caputo curve in the y-direction
# Values monotonically decrease, as in with distance from fluvial source
# Modified from what is in fluvialProfileModule.R to always put the left coordinate before the right coordinate (as in xRescale)
# Used by generateInitialProfile() and profileForShore(), ultimately by fillBasin()
yRescale <- function(y, leftY, rightY) {
	maxY <- max(y)
	minY <- min(y)
	rescaled <- y * (leftY - rightY) / (maxY - minY) + rightY
	rescaled
}

# Calculate the integrated sedimentation over the profile (in m*km), given initial and subsequent topographic profiles
# Used by profileForShore(), in turn by optimalProfile() and fillBasin()
sedimentVolume <- function(initial, subsequent, deltaX, includeErosion=TRUE) {
	sediment <- subsequent$y - initial$y
	if (includeErosion == FALSE) {
		sediment[sediment < 0] <- 0
	}
	totalSediment <- sum(sediment, na.rm=TRUE) * deltaX   # units will be m*km
	totalSediment
}

# Generate the starting nonmarine and marine profiles for a basin simulation
# Used by fillBasin()
generateInitialProfile <- function(fallLineX=0, fallLineY=150, shoreX=200, shoreY=0, deltaWidth=100, marginWidth=400, deltaToeY=-100, nonMarAlpha=0.5, marineAlpha=2.0, deltaX=1) {
	# Lateral positions and lengths are in km, including fallLineX, shoreX, deltaWidth,
	#   marginWidth, and deltaX. deltaX is the spatial resolution of the model.
	# Vertical positions are elevations above sea level, in m, including fallLineY, shoreY, and deltaToeY
	# nonMarAlpha and marineAlpha are used to set the concavity of the Caputo profiles

	# Caputo profiles originally scaled from 0 to 1
	#  rescaled to the length of the x vector (determined by marginWidth and deltaX)
	#  interpolated at end to produce points at the fixed intervals of deltaX km
	x <- seq(0, 1, 0.001)      # to make finely resolved Caputo profiles
	# Changing 0.001 to 0.0001 increases the run time by 10x

	# nonmarine Caputo profile
	h <- height(x, alpha=nonMarAlpha)
	topoNonmarineX <- xRescale(x, leftX=fallLineX, rightX=shoreX)
	topoNonmarineY <- yRescale(h, leftY=fallLineY, rightY=shoreY)

	# delta Caputo profile
	h <- height(x, alpha=marineAlpha)
	topoDeltaX <- xRescale(x, leftX=shoreX, rightX=shoreX + deltaWidth)
	topoDeltaY <- yRescale(h, leftY=shoreY, rightY=deltaToeY)

	# pro-deltaic shelf (non Caputo; flat at start of simulation)
	topoShelfX <- seq(shoreX + deltaWidth, marginWidth, 1)
	topoShelfY <- rep(deltaToeY, length(topoShelfX))

	# combine the three profile segments, removing duplicated points where profiles meet
	topoX <- c(topoNonmarineX, topoDeltaX[-1], topoShelfX[-1])
	topoY <- c(topoNonmarineY, topoDeltaY[-1], topoShelfY[-1])

	# calculate interpolated profile to make the points fall on the x-grid of the model
	interpolated <- stats::approx(topoX, topoY, xout=seq(fallLineX, marginWidth, deltaX))
	return(interpolated)
}

# Calculates the profile for a specified shore position. Code should be identical to profileFit, except that it returns the interpolated profile
# Used by optimalProfile()
profileForShore <- function(deflected, shoreX, deltaWidth, marginWidth, deltaX, nonMarAlpha, marineAlpha, targetSedVolume, fallLineY) {
	fallLineX <- deflected$x[1]
	deltaToeX <- shoreX + deltaWidth
	if (deltaToeX > marginWidth) {
		warning("Delta toe extends beyond the margin. Increase marginWidth.", call.=FALSE, immediate.=TRUE)
		print(paste('deltaToeX: ', deltaToeX, '     marginWidth: ', marginWidth, '\n'))
	}

	# y coordinates are elevations in m
	# fallLineY <- deflected$y[1]   TODO: Delete this line
	shoreY <- 0
	deltaToeY <- deflected$y[indexForPosition(positionKm=deltaToeX, marginWidth=marginWidth, deltaX=deltaX)]

	# Caputo profiles originally scaled from 0 to 1
	#  rescaled to the length of the x vector (determined by marginWidth and deltaX)
	#  interpolated at end to produce points at the fixed intervals of deltaX km
	x <- seq(0, 1, 0.001)     # to make finely resolved Caputo profiles
	# Changing 0.001 to 0.0001 increases the run time by 10x

	# nonmarine Caputo profile
	h <- height(x, alpha=nonMarAlpha)
	topoNonmarineX <- xRescale(x, leftX=fallLineX, rightX=shoreX)
	topoNonmarineY <- yRescale(h, leftY=fallLineY, rightY=shoreY)

	# delta Caputo profile
	h <- height(x, alpha=marineAlpha)
	topoDeltaX <- xRescale(x, leftX=shoreX, rightX=deltaToeX)
	topoDeltaY <- yRescale(h, leftY=shoreY, rightY=deltaToeY)

	# pro-deltaic shelf (non Caputo; only remnant topography)
	topoShelfX <- seq(ceiling(deltaToeX), marginWidth, deltaX)
	deltaToeIndex <- indexForPosition(positionKm=ceiling(deltaToeX), marginWidth=marginWidth, deltaX=deltaX)
	marginWidthIndex <- indexForPosition(positionKm=marginWidth, marginWidth=marginWidth, deltaX=deltaX)
	topoShelfY <- deflected$y[deltaToeIndex:marginWidthIndex]

	# combine the three profile segments, removing duplicated points where nonmarine and delta profiles meet
	topoX <- c(topoNonmarineX, topoDeltaX[-1], topoShelfX)
	topoY <- c(topoNonmarineY, topoDeltaY[-1], topoShelfY)

	# check for NAs
	if (any(is.na(topoX)) == TRUE) {
		warning("Unexpected NA in topoX.", call.=FALSE, immediate.=TRUE)
		# This will cause a problem in the stats::approx() function
	}

	if (any(is.na(topoY)) == TRUE) {
		warning("Unexpected NA in topoX.", call.=FALSE, immediate.=TRUE)
		# These points will be ignored in the stats::approx() function, but would be good to know if this happens.
	}

	# calculate interpolated profile to make the points fall on the x-grid of the model
	interpPoints <- seq(fallLineX, marginWidth, deltaX)
	interpolated <- stats::approx(topoX, topoY, xout=interpPoints, rule=2, ties='ordered')

	# prevent marine erosion
	if (any(is.na(interpolated$y)) == TRUE) {
		warning("Unexpected NA in interpolation. profileForShore()", call.=FALSE, immediate.=TRUE)
		# Originally thought to be attempting to erode through the floor of the shelf, thought this was caused by too high of a magnitude of sea-level fall. Now known to occur in a subsidence-only run. It's unclear why, but the problem is at the highest (rightmost X-point), where the interpolated y value becomes NA. The fix below is simply to replace that NA with the topoY value at that point (they ought to match). It would be better if I understood the root cause of this NA value, and it's not clear from the documentation.
		# save(fallLineX, marginWidth, deltaX, interpPoints, topoX, topoY, interpolated, file='~/Desktop/interpolatedDebug.RData')
		# attempt to fix interpolation error on final point
		if (is.na(interpolated$y[length(interpolated$y)]) == TRUE) {
			interpolated$y[length(interpolated$y)] = topoY[length(topoY)]
		}
	}
	# prevent marine erosion
	marine <- interpolated$y < 0
	if (length(interpolated$y) != length(deflected$y)) {
		warning("Length of interpolated$y and deflected$y are unequal in profileForShore().", call.=FALSE, immediate.=TRUE)
		print(paste('interpolated$y:', length(interpolated$y), ', deflected$y:', length(deflected$y)));
	}
	erosion <- interpolated$y < deflected$y
	interpolated$y[marine & erosion] <- deflected$y[marine & erosion]

	sedVolume <- sedimentVolume(initial=deflected, subsequent=interpolated, deltaX=deltaX, includeErosion=FALSE)
	volumeDeviation <- abs(sedVolume - targetSedVolume)

	results <- list(topoProfile=interpolated, sedVolume=sedVolume, fit=volumeDeviation)
	return(results)
}

# Calculates how well the profile satisfies the targetSedVolume, by returning the absolute value of the difference between the actual sediment volume and the target sediment volume.
# Used by optimalProfile(), which optimizes this function
profileFit <- function(deflected, shoreX, deltaWidth, marginWidth, deltaX, nonMarAlpha, marineAlpha, targetSedVolume, fallLineY) {
	topoProfile <- profileForShore(deflected=deflected, shoreX=shoreX, deltaWidth=deltaWidth, marginWidth=marginWidth, deltaX=deltaX, nonMarAlpha=nonMarAlpha, marineAlpha=marineAlpha, targetSedVolume=targetSedVolume, fallLineY=fallLineY)
	return(topoProfile$fit)
}

# Calculates an optimal profile, where the shore is placed to most closely achieve the target sediment volume
# Used by fillBasin()
optimalProfile <- function(deflected, previousShoreX, deltaWidth, marginWidth, deltaX, searchWidth, targetSedVolume, nonMarAlpha, marineAlpha, fallLineY) {
	# deflected is the previous land surface, adjusted for eustasy and subsidence
	# previousShoreX, deltaWidth, marginWidth, deltaX, searchWidth all in km
	# targetSedVolume in m * km
	# nonMarAlpha and marineAlpha are Caputo shape parameters

	# initial pass to find approximate shoreline, to nearest km
	shore <- seq(round(previousShoreX) - searchWidth, round(previousShoreX) + searchWidth, 1.0)
	sedVolumes <- rep(0, length(shore))
	for (i in 1:length(shore)) {
		theProfile <- profileForShore(deflected=deflected, shoreX=shore[i], deltaWidth=deltaWidth, marginWidth=marginWidth, deltaX=deltaX, nonMarAlpha=nonMarAlpha, marineAlpha=marineAlpha, targetSedVolume=targetSedVolume, fallLineY=fallLineY)
		sedVolumes[i] <- theProfile$sedVolume
	}

	volumeDeviation <- abs(sedVolumes - targetSedVolume)
	bestShoreIndex <- which(volumeDeviation == min(volumeDeviation))
	bestShoreIndex <- bestShoreWarnings(bestShoreIndex=bestShoreIndex, sedVolumes=sedVolumes, shore=shore)
	bestShore <- shore[bestShoreIndex]
	windowPadding <- 1.0   # was 2.0 (if this is too large, local minima may result)
	searchWindow <- c(bestShore-windowPadding, bestShore+windowPadding)

	# use initial pass to search for optimal shore
	optimalFit <- stats::optimize(profileFit, searchWindow, deflected=deflected, deltaWidth=deltaWidth, marginWidth=marginWidth, deltaX=deltaX, nonMarAlpha=nonMarAlpha, marineAlpha=marineAlpha, targetSedVolume=targetSedVolume, fallLineY=fallLineY)
	optimalShore <- optimalFit$minimum
	optimizedProfile <- profileForShore(deflected, optimalShore, deltaWidth, marginWidth, deltaX, nonMarAlpha, marineAlpha, targetSedVolume, fallLineY)
	optimalTopo <- optimizedProfile$topoProfile
	optimalVolume <- optimizedProfile$sedVolume

	# in case a local minimum was found, decrease the search window and try again
	while (optimalVolume / targetSedVolume < 0.5) {
		message("local minimum encountered; attempting to correct")
		windowPadding <- windowPadding / 2
		searchWindow <- c(bestShore-windowPadding, bestShore+windowPadding)
		optimalFit <- stats::optimize(profileFit, searchWindow, deflected=deflected, deltaWidth=deltaWidth, marginWidth=marginWidth, deltaX=deltaX, nonMarAlpha=nonMarAlpha, marineAlpha=marineAlpha, targetSedVolume=targetSedVolume, fallLineY=fallLineY)
		optimalShore <- optimalFit$minimum
		optimizedProfile <- profileForShore(deflected, optimalShore, deltaWidth, marginWidth, deltaX, nonMarAlpha, marineAlpha, targetSedVolume, fallLineY)
		optimalTopo <- optimizedProfile$topoProfile
		optimalVolume <- optimizedProfile$sedVolume

		if (windowPadding < 0.1) {
		    stop("ERROR: persistent local minimum encountered; unable to correct")
		}
	}

	# if a local minimum persists, display a message
	if (optimalVolume / targetSedVolume < 0.5) {
		message("WARNING: persistent local minimum encountered; unable to correct")
	}

	# DEBUGGING ONLY: IDENTIFY PARTITIONING SPIKES / NEAR-ZERO SEDIMENTATION
	# REFLECTS stats::optimize() FINDING A LOCAL MINIMUM
#	if (optimalVolume / targetSedVolume < 0.5) {
#		message(paste("targetSedVolume:", targetSedVolume, "; optimalVolume:", optimalVolume))
#		fileName = paste(round(stats::rnorm(1)*100+1000), ".RData", sep="")
#		save(profileFit, searchWindow, deflected, deltaWidth, marginWidth, deltaX, nonMarAlpha, marineAlpha, targetSedVolume, fallLineY=fallLineY, shore, sedVolumes, volumeDeviation, previousShoreX, bestShoreIndex, bestShore, windowPadding, optimalFit, optimalShore, optimizedProfile, optimalTopo, optimalVolume, file=fileName)
#	}

	results <- list(shore=optimalShore, sedimentVolume=optimalVolume, fittedProfile=optimalTopo)
	return(results)
}




# ---------------------------------------------------------------------------------------
# ERROR CHECKING

# Called by fillBasin()
# Verifies that all objects have the same number of time steps and spatial points
integrityCheck <- function(geometry, subsidence, eustasy, sediment) {
	haltSimulation <- FALSE

	# Time step check
	geometryTimeSteps <- geometry$duration / geometry$timeStep + 1
	eustasyTimeSteps <- length(eustasy$timeSeries$timePoint)
	subsidenceTimeSteps <- nrow(subsidence$rates)
	sedimentTimeSteps <- length(sediment$timeSeries$timePoint)
	if (eustasyTimeSteps != geometryTimeSteps) {
		message(paste("Error: Number of time steps in eustasy (", eustasyTimeSteps, ") does not match the number of time steps implied by geometry (", geometryTimeSteps, ")\r", sep=""))
		haltSimulation <- TRUE
	}
	if (subsidenceTimeSteps != geometryTimeSteps) {
		message(paste("Error: Number of time steps in subsidence (", subsidenceTimeSteps, ") does not match the number of time steps implied by geometry (", geometryTimeSteps, ")\r", sep=""))
		haltSimulation <- TRUE
	}
	if (sedimentTimeSteps != geometryTimeSteps) {
		message(paste("Error: Number of time steps in sediment (", sedimentTimeSteps, ") does not match the number of time steps implied by geometry (", geometryTimeSteps, ")\r", sep=""))
		haltSimulation <- TRUE
	}

	# Spatial points check
	geometryPoints <- (geometry$marginWidth - geometry$fallLineX) / geometry$deltaX + 1
	subsidencePoints <- ncol(subsidence$rates)
	if (subsidencePoints != geometryPoints) {
		message(paste("Error: Number of spatial points for subsidence (", subsidencePoints, ") does not match the number of spatial points implied by geometry (", geometryPoints, ")\r", sep=""))
		haltSimulation <- TRUE
	}

	if (haltSimulation == TRUE) {
		stop("Model cannot be run until these issues are fixed.\n", call.=FALSE)
	}

}

# Various tests for the best shore, called by optimalProfile() to simplify its code. These issues may have been solved, as I haven't seen these warnings in a while.
# Used by optimalProfile()
bestShoreWarnings <- function(bestShoreIndex, sedVolumes, shore) {
	if (length(bestShoreIndex) > 1) {
		if (stats::median(sedVolumes[bestShoreIndex]) < 0.0001) {
			# In some cases, multiple best shore indices are found, all with a sediment
			# volume of 0. If the target sediment volume is small, the next highest index
			# might have too large of a sediment volume. Rather than deposit no sediment,
			# the next-highest index will be selected as the best shore position
			if (max(bestShoreIndex) < length(shore)) {
				bestShoreIndex <- max(bestShoreIndex) + 1
				# warning is turned off in version 0.4.0 because it rarely indicated an actual problem
				# warning("bestShoreIndex too long - case 1; usually not a problem.", call.=FALSE, immediate.=FALSE) # rather common, but not usually a problem
			} else {
				bestShoreIndex <- max(bestShoreIndex)
				warning("Unable to find best shore, results likely wrong. Decrease deltaX, and if necessary, increase searchWidth.", call.=FALSE, immediate.=TRUE)
			}
		} else {
			# For other ties, choose (arbitrarily) the left-most of these
			bestShoreIndex <- bestShoreIndex[1]
			warning("bestShoreIndex too long - case 2; usually not a problem.", call.=FALSE, immediate.=TRUE)
		}
	} else if (length(bestShoreIndex) < 1) {
		warning("No best shore was found", call.=FALSE, immediate.=TRUE)
	}

	if (shore[1] <= 0) {
		warning("Shore is too close to the left edge of the basin; run is likely unusable.", call.=FALSE, immediate.=TRUE)
	} else if (bestShoreIndex == 1 | bestShoreIndex == length(shore)) {
		warning("Best shore is at edge of search window; searchWidth is too small.", call.=FALSE, immediate.=TRUE)
	}

	return(bestShoreIndex)
}
