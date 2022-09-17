# Occurrences Module

# ---------------------------------------------------------------------------------------
## SIMULATE FOSSIL OCCURRENCES

# Generate occurrences from stratigraphic history and species data
generateOccurrences <- function(stratColumn, species, sampleSpacing=0.5) {
	attach(stratColumn)  # modelTime, elevation, stratPosition, facies
	                                  # unused: thickness
	attach(species)  # id, ancestor, origination, extinction, PE, ET, PA, Aff
	
	numSpecies <- length(id)
	
	# vectors to hold results; oversize, then cull at end if necessary
	# age and position are equivalent to modelTime and elevation, will be renamed when
	#   data frame is created, but not here to avoid collision with attached strat data
	maxLength <- 100000
	storage <- vector(mode='numeric', length=maxLength)
	age       <- storage
	position  <- storage
	speciesId <- storage
	occElevation <- storage
	occFacies    <- storage
	
	collectionHorizons <- seq(0, max(stratPosition), sampleSpacing)
	occurrenceNum <- 1
	
	for (horizon in collectionHorizons) {
		horizonTime <- modelTimeForHorizon(stratColumn, horizon, oldestAge=FALSE)
		horizonElevation <- elevation[horizonTime == modelTime]
		horizonFacies <- facies[horizonTime == modelTime]
		for (i in 1:numSpecies) {
			if (horizonTime >= origination[i] && horizonTime <= extinction[i]) { # extant
				multiplier <- preservation(affinity=Aff[i], facies=horizonFacies)
				probCollection <- multiplier * probCollection(PE=PE[i], ET=ET[i], PA=PA[i], elevation=horizonElevation)						
				random <- runif(1)
				if (random <= probCollection) {	                         # found
					age[occurrenceNum] <- horizonTime
					position[occurrenceNum] <- horizon
					speciesId[occurrenceNum] <- id[i]
					occElevation[occurrenceNum] <- horizonElevation
					occFacies[occurrenceNum] <- horizonFacies
					occurrenceNum <- occurrenceNum + 1
					if (occurrenceNum%%maxLength == 0) {
						age <- c(age, storage)
						position <- c(position, storage)
						speciesId <- c(speciesId, storage)
						occElevation<- c(horizonElevation, storage)
						occFacies <- c(horizonFacies, storage)
						message(paste("Resizing storage to hold up to", length(speciesId), "occurrences. If this happens multiple times, this may be a sign that the run needs to be terminated."))
					}
				}
			}
		}
	}	
	
	detach(stratColumn)
	detach(species)
	occurrences <- data.frame(modelTime=age, stratPosition=position, speciesId, elevation=occElevation, facies=occFacies)
	# trim any rows with no data
	occurrences <- occurrences[-which(occurrences$speciesId == 0), ]
	message(paste("Generated", nrow(occurrences), "occurrences"))
	occurrences
}

# Generate occurrences, but by stepping through time and finding the next horizon above
# the previous horizon. If sampleSpacing is less than channelDepth, channels will be 
# sampled only once by this method, causing denser sampling in floodplain facies than
# channels, potentially useful if you think of channels only containing fossils in their
# basal lag.
# [TEST]
generateOccurrencesByTime <- function(stratColumn, species, sampleSpacing=0.5) {
	attach(stratColumn)  # modelTime, elevation, stratPosition
	                                  # unused: thickness, facies
	attach(species)  # id, ancestor, origination, extinction, PE, ET, PA
	
	numSpecies <- length(id)
	
	# vectors to hold results; oversize, then cull at end if necessary
	# age and position are equivalent to modelTime and elevation, will be renamed when
	#   data frame is created, but not here to avoid collision with attached strat data
	maxLength <- 100000
	storage <- vector(mode='numeric', length=maxLength)
	age       <- storage
	position  <- storage
	speciesId <- storage
	
	occurrenceNum <- 1
	previousSamplePosition <- min(stratPosition) - 2*sampleSpacing
		# insures lowest horizon will be used
	for (timeStep in 1:length(modelTime)) {
		if (stratPosition[timeStep] >= previousSamplePosition + sampleSpacing) { # get sample
			previousSamplePosition <- stratPosition[timeStep]
			for (i in 1:numSpecies) {
				if (modelTime[timeStep] >= origination[i] && modelTime[timeStep] <= extinction[i]) { # extant
					probCollection <- probCollection(PE=PE[i], ET=ET[i], PA=PA[i], elevation=elevation[timeStep])						
					random <- runif(1)
					if (random <= probCollection) {	                         # found
						age[occurrenceNum] <- modelTime[timeStep]
						position[occurrenceNum] <- stratPosition[timeStep]
						speciesId[occurrenceNum] <- id[i]
						occurrenceNum <- occurrenceNum + 1
						if (occurrenceNum%%maxLength == 0) {
							age <- c(age, storage)
							position <- c(position, storage)
							speciesId <- c(speciesId, storage)
							message(paste("Resizing storage to hold up to", length(speciesId), "occurrences. If this happens multiple times, this may be a sign that the run needs to be terminated."))

						}
					}
				}
			}
		}
		
	}
	
	detach(stratColumn)
	detach(species)
	occurrences <- as.data.frame(cbind(modelTime=age, stratPosition=position, speciesId))
	# trim any rows with no data
	occurrences <- occurrences[-which(occurrences$speciesId == 0), ]
	message(paste("Generated", nrow(occurrences), "occurrences"))
	occurrences
}

# ---------------------------------------------------------------------------------------
# CONVERSION FROM AFFINITY TO MULTIPLIER

preservation <- function(affinity, facies) {
	multiplier <- 1
	if (facies == 'channel') {
		if (affinity > 0) {  # species found preferentially in floodplains
			multiplier <- 1 - affinity
		}
	} else if (facies == 'floodplain') {
		if (affinity < 0) {  # species found preferentially in channels
			multiplier <- 1 + affinity
		}
	} else {   # unspecified, marine, paleosol, and unconformity
		multiplier <- 0
	}
	return(multiplier)
}



# ---------------------------------------------------------------------------------------
## ANALYSIS OF FOSSIL OCCURRENCES

# Find the first and last occurrence of every species in an occurrence file
# Used in occurrencePlot
fadAndLad <- function(occurrences) {
	# occurrences should have 3 columns: modelTime, stratPosition, speciesId
	occurringSpecies <- sort(unique(occurrences$speciesId))
	fad <- vector(mode='numeric', length=length(occurringSpecies))
	lad <- vector(mode='numeric', length=length(occurringSpecies))
	fadTime <- vector(mode='numeric', length=length(occurringSpecies))
	ladTime <- vector(mode='numeric', length=length(occurringSpecies))
	numOccurrences <- vector(mode='numeric', length=length(occurringSpecies))
	
	for (i in 1:length(occurringSpecies)) {
		speciesOccurrence <- occurrences[occurrences$speciesId == occurringSpecies[i], ]
		fad[i] <- min(speciesOccurrence$stratPosition)
		lad[i] <- max(speciesOccurrence$stratPosition)
		fadTime[i] <- min(speciesOccurrence$modelTime)
		ladTime[i] <- max(speciesOccurrence$modelTime)
		numOccurrences[i] <- length(speciesOccurrence$stratPosition)
	}
	
	fadAndLad <- as.data.frame(cbind(occurringSpecies, fad, lad, fadTime, ladTime, numOccurrences))
	return(fadAndLad)
}

# Calculate the percent of species that leave a record, are singletons, and are two-timers
completeness <- function(occurrences, species) {
	totalSpecies <- nrow(species)	
	tallies <- table(occurrences$speciesId)
	occurringSpecies <- length(tallies)
	singletons <- length(tallies[tallies == 1])
	twoTimers  <- length(tallies[tallies == 2])
	
	percentOccurrence <- round(occurringSpecies/totalSpecies*100, 1)
	percentSingletons <- round(singletons/totalSpecies*100, 1)
	percentTwoTimers  <- round(twoTimers/totalSpecies*100, 1)
	
	message(paste("There are ", totalSpecies, " species in the simulation, of which", sep=''))
	message(paste("  ", occurringSpecies, " (", percentOccurrence, "%) left a fossil record, ", sep=''))
	message(paste("  ", singletons, " (", percentSingletons, "%) occur only once, and", sep=''))
	message(paste("  ", twoTimers, " (", percentTwoTimers, "%) occur only twice.", sep=''))
	message(" ")
	
	results <- list(totalSpecies=totalSpecies, occurringSpecies=occurringSpecies, singletons=singletons, twoTimers=twoTimers)
	return(results)
}

# Calculate the range offset for fads and lads (difference in age between origination and fad, or between extinction and lad)
rangeOffset <- function(occurrences, species) {
	occurringSpecies <- sort(unique(occurrences$speciesId))
	numOccurringSpecies <- length(occurringSpecies)
	fad <- vector(mode='numeric', length=numOccurringSpecies)
	lad <- vector(mode='numeric', length=numOccurringSpecies)
	fadOffset <- vector(mode='numeric', length=numOccurringSpecies)
	ladOffset <- vector(mode='numeric', length=numOccurringSpecies)
	
	for (i in 1:numOccurringSpecies) {
		currentSpeciesId <- occurringSpecies[i]
		speciesOccurrence <- occurrences[occurrences$speciesId == currentSpeciesId, ]
		fad[i] <- min(speciesOccurrence$stratPosition)
		lad[i] <- max(speciesOccurrence$stratPosition)
		
		fadTime <- min(speciesOccurrence$modelTime)
		origination <- species$origination[currentSpeciesId]
		fadOffset[i] <- fadTime - origination
		
		ladTime <- max(speciesOccurrence$modelTime)
		extinction <- species$extinction[currentSpeciesId]
		ladOffset[i] <- extinction - ladTime
	}
	
	results <- data.frame(occurringSpecies, fad, lad, fadOffset, ladOffset)
	return(results)
}


# ---------------------------------------------------------------------------------------
## OCCURRENCE PLOTS

# Plot a range chart with occurrences
occurrencePlot <- function(occurrences, stratColumn, species, peakExtTimeMy=NA, extDurationMy=NA, extinctionColor='red', orderedByLad=TRUE, yAxisLabels=TRUE) {
	# determine fads and lads
	fadsLads <- fadAndLad(occurrences)
	if (orderedByLad == TRUE) {
		fadsLads <- fadsLads[c(order(fadsLads$lad)),]
	} else {
		fadsLads <- fadsLads[c(order(fadsLads$fad)),]
	}
	rows <- seq(1:length(fadsLads$fad))
	
	occurrenceColor <- rgb(0, 0, 1, 0.5)
	ylab = 'Stratigraphic Position (m)'
	if (yAxisLabels == FALSE) {
		ylab = ''
	}
	
	plot(rows, fadsLads$fad, ylim=c(min(stratColumn$stratPosition), max(stratColumn$stratPosition)), xlab='Species', ylab=ylab, las=1, type='n', frame.plot=yAxisLabels, axes=yAxisLabels)
	if (yAxisLabels == FALSE) {
		axis(1)
	}
	
	# if there is a mass extinction, add a box for the interval and a line for the peak
	if (massExtinctionOccurred(peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)) {
		addExtinctionBox(stratColumn=stratColumn, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy, leftEdge=min(rows), rightEdge=max(rows), extinctionColor=extinctionColor)
	}
	
	# times of origination and extinction
	for (i in 1:length(fadsLads$occurringSpecies)) {
		theSpecies <- fadsLads$occurringSpecies[i]
		timeOrigination <- species$origination[theSpecies]
		timeExtinction  <- species$extinction[theSpecies]
		horizonOrigination <- stratPositionForAge(timeOrigination, stratColumn)
		horizonExtinction <- stratPositionForAge(timeExtinction, stratColumn)
		points(i, horizonOrigination, pch=3, cex=0.5, col='gray')
		points(i, horizonExtinction, pch=3, cex=0.5, col='gray')
		segments(i, horizonOrigination, i, horizonExtinction, lwd=0.5, lty='dotted', col='gray')
	}

	# add ranges
	segments(rows, fadsLads$fad, rows, fadsLads$lad, col='blue', lwd=1)

	# add the occurrences to the plot
	for (i in 1:length(fadsLads$occurringSpecies)) {
		theSpecies <- fadsLads$occurringSpecies[i]
		horizons <- occurrences$stratPosition[occurrences$speciesId == theSpecies]
		x <- rep(i, length(horizons))
		points(x, horizons, col=occurrenceColor, pch=16, cex=0.7)
	}

	# singletons
	points(rows[fadsLads$numOccurrences == 1], fadsLads$fad[fadsLads$numOccurrences == 1], pch=16, cex=0.7, col='white')
	points(rows[fadsLads$numOccurrences == 1], fadsLads$fad[fadsLads$numOccurrences == 1], pch=1, cex=0.8, col=occurrenceColor)
}

# Plot all occurrences, along with a stratigraphic column and the history of elevation
occurrenceColumnElevationPlot <- function(occurrences, stratColumn, species, peakExtTimeMy=NA, extDurationMy=NA, orderedByLad=TRUE) {
	opar <- par(no.readonly = TRUE)
	leftDivider <- 0.20   # separates ranges from elevation plot
	rightDivider <- 0.80  # separates elevation plot from strat column
	
	# Left panel: stratigraphic column
	par(fig=c(0, leftDivider, 0, 1), mar=c(5, 4, 4, 0) + 0.1)
	stratColumnPlot(stratColumn)
	
	# Middle panel: Range chart
	par(fig=c(leftDivider, rightDivider, 0, 1), mar=c(5, 0, 4, 0) + 0.1, new=TRUE)
	occurrencePlot(occurrences=occurrences, stratColumn=stratColumn, species=species, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy, orderedByLad=orderedByLad, yAxisLabels=FALSE)
	
	# Right panel: elevation plot
	par(fig=c(rightDivider, 1, 0, 1), mar=c(5, 0, 4, 2) + 0.1, new=TRUE)
	elevationHistoryPlot(stratColumn, yAxisLabels=FALSE)
	
	par(opar)
}

# Plot bar chart of fads or lads in each horizon from occurrence data
fadLadPlot <- function(type=c('fad', 'lad'), occurrences, stratColumn, sampleSpacing=0.5, xMax=NA, peakExtTimeMy=NA, extDurationMy=NA, extinctionColor='red') {
	fadsLads <- fadAndLad(occurrences)
	bins <- 0
	fadLadHist <- 0
	xlab <- ''
	
	if (type == 'fad') {
		upperLimit <- ceiling(max(fadsLads$fad))
		if (upperLimit %% sampleSpacing > 0) upperLimit <- upperLimit + sampleSpacing
		lowerLimit <- floor(min(fadsLads$fad))
		if (lowerLimit %% sampleSpacing > 0) lowerLimit <- 0 # base of section
		bins <- seq(lowerLimit, upperLimit, sampleSpacing)
		fadLadHist <- hist(fadsLads$fad, breaks=bins, plot=FALSE)
		xlab <- 'First occurrences'
	} else if (type == 'lad') {
		upperLimit <- ceiling(max(fadsLads$lad))
		if (upperLimit %% sampleSpacing > 0) upperLimit <- upperLimit + sampleSpacing
		lowerLimit <- floor(min(fadsLads$lad))
		if (lowerLimit %% sampleSpacing > 0) lowerLimit <- 0 # base of section
		bins <- seq(lowerLimit, upperLimit, sampleSpacing)
		fadLadHist <- hist(fadsLads$lad, breaks=bins, plot=FALSE)
		xlab <- 'Last occurrences'
	} else {
		stop('unrecognized plot type')
	}
	
	barHalfWidth <- sampleSpacing/2
	if (!is.na(xMax) &  xMax < max(fadLadHist$counts)) {
		fadLadHist$counts[fadLadHist$counts > xMax] <- xMax
	}
	xlim = range(fadLadHist$counts)
	if (!is.na(xMax)) {
		xlim = c(0, xMax)
	}
	ylim = c(min(stratColumn$stratPosition), max(stratColumn$stratPosition))
	plot(fadLadHist$counts, fadLadHist$mids, type='n', xlim=xlim, ylim=ylim, xlab=xlab, ylab='Stratigraphic Position (m)', las=1)
	
	if (massExtinctionOccurred(peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)) {
		addExtinctionBox(stratColumn=stratColumn, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy, leftEdge=min(fadLadHist$counts), rightEdge=max(fadLadHist$counts), extinctionColor=extinctionColor)
	}
	
	nonZero <- fadLadHist$counts > 0
	rect(0, fadLadHist$mids[nonZero]-barHalfWidth, fadLadHist$counts[nonZero], fadLadHist$mids[nonZero] + barHalfWidth, col='gray')
}

# Combined plot showing number of fads or lads per horizon, elevation history, and sedimentation rate history
fadLadElevationSedPlot <- function(type=c('fad', 'lad'), occurrences, stratColumn, sampleSpacing=0.5, stratBin=10, peakExtTimeMy=NA, extDurationMy=NA) {
	opar <- par(no.readonly = TRUE)
	par(cex.axis=0.8)
	leftDivider <- 0.395  # picked by trial and error to give correct proportions
	rightDivider <- 0.675
	
	# left panel: first occurrences
	par(fig=c(0, leftDivider, 0, 1), mar=c(5, 4, 4, 0.5) + 0.1)
	fadLadPlot(type=type, occurrences=occurrences, stratColumn=stratColumn, sampleSpacing=sampleSpacing, xMax=NA, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)
	
	# middle panel: elevation
	par(fig=c(leftDivider, rightDivider, 0, 1), mar=c(5, 0, 4, 0.5) + 0.1, new=TRUE)
	elevationHistoryPlot(stratColumn, yAxisLabels=FALSE)
	
	# right panel: sedimentation rate
	par(fig=c(rightDivider, 1, 0, 1), mar=c(5, 0, 4, 2) + 0.1, new=TRUE)
	sedRateHistoryPlot(stratColumn, stratBin=stratBin, yAxisLabels= FALSE)
	
	par(opar)
}


# Plot diversity vs. elevation, used as a diagnostic plot to see if diversity covaries with elevation (it generally shouldn't)
# In this function, a species is considered to occur if its probability of collection exceeds a specified threshold (pCrit)
# elevationDiversityResamplePlot() is the preferred way to do this, as it adds a confidence interval, but it is also slower
elevationDiversityPlot <- function(species, minElevation=0, maxElevation=150, pCrit=0.1) {
	attach(species)
	numSpecies <- nrow(species)
	
	elevation <- seq(minElevation, maxElevation, 1)
	diversity <- rep(0, length(elevation))
	
	for (elev in 1:length(elevation)) {
		elevDiversity <- 0
		for (i in 1:numSpecies) {
			probCollection <- probCollection(PE=PE[i], ET=ET[i], PA=PA[i], elevation=elevation[elev])
			if (probCollection >= pCrit) {	                         # found
				elevDiversity <- elevDiversity + 1
			}
		}
		diversity[elev] <- elevDiversity
	}
	
	detach(species)
	plot(elevation, diversity, ylim=c(0, max(diversity)), xlab='Elevation (m)', ylab='Richness (number of species)', type='l', lwd=2, las=1)
}

# Diagnostic plot to evaluate whether diversity changes with elevation (it generally shouldn't)
# In this function, a species is considered to occur if a random number is less than its probability of collection.
# This is done over a number of trials to place a confidence interval on the observed diversity at any elevation.
# Function is slow, owing to the three nested loops.
elevationDiversityResamplePlot <- function(species, minElevation=0, maxElevation=150, time=0, numTrials=100) {
	selectedSpecies <- species[which(species$origination <= time & species$extinction >= time), ]
	numSpecies <- nrow(selectedSpecies)
	attach(selectedSpecies)
	
	elevation <- seq(minElevation, maxElevation, 1)
	diversity50 <- rep(0, length(elevation)) # mean
	diversity25 <- rep(0, length(elevation)) # 25% quantile
	diversity75 <- rep(0, length(elevation)) # 75% quantile
	
	for (elev in 1:length(elevation)) {
		trialDiversity <- rep(0, numTrials)
		for (trial in 1:numTrials) {
			iterationDiversity <- 0
			for (i in 1:numSpecies) {
				random <- runif(1)
				probCollection <- probCollection(PE=PE[i], ET=ET[i], PA=PA[i], elevation=elevation[elev])
				if (random <= probCollection) {	                         # found
					iterationDiversity <- iterationDiversity + 1
				}
			}
			trialDiversity[trial] <- iterationDiversity
		}
		diversity50[elev] <- quantile(trialDiversity, 0.50)
		diversity25[elev] <- quantile(trialDiversity, 0.25)
		diversity75[elev] <- quantile(trialDiversity, 0.75)
	}

	detach(selectedSpecies)
	yLimits <- c(0, max(diversity75))
	plot(elevation, diversity50, ylim=yLimits, xlab='Elevation (m)', ylab='Richness (number of species)', type='n', lwd=2, las=1)

	polyX <- c(elevation, rev(elevation))
	polyY <- c(diversity25, rev(diversity75))
	polygon(polyX, polyY, col='lightgray', border='lightgray')
	
	points(elevation, diversity50, type='l', lwd=2, col='black')
}




# ---------------------------------------------------------------------------------------
## UTILITY FUNCTIONS

# Probability of collection for a species at a given elevation
# Used by generateOccurrences(), generateOccurrencesByTime(), and the elevation-diversity plots
probCollection <- function(PE, ET, PA, elevation) {
	# PA is on a percent scale: 0-100
	probCollection <- (PA * exp(-((elevation - PE)^2)/(2 * ET^2))) / 100
	if (elevation <= 0) {           # nonmarine (marine would have a negative elevation)
		probCollection <- 0
	}
	probCollection
}

# Find horizons corresponding to the onset, peak, and end of a mass extinction
# Used by addExtinctionBox()
extinctionHorizons <- function(stratColumn, peakExtTimeMy, extDurationMy) {
	extStart <- peakExtTimeMy - extDurationMy/2
	extEnd <- peakExtTimeMy + extDurationMy/2
	
	timeDiff <- stratColumn$modelTime - extEnd
	index1 <- which(timeDiff == min(abs(timeDiff)))
	index2 <- which(timeDiff == -min(abs(timeDiff)))
	endIndex <- max(index1, index2)
	extEndHorizon <- stratColumn$stratPosition[endIndex]

	timeDiff <- stratColumn$modelTime - extStart
	index1 <- which(timeDiff == min(abs(timeDiff)))
	index2 <- which(timeDiff == -min(abs(timeDiff)))
	startIndex <- max(index1, index2)
	extStartHorizon <- stratColumn$stratPosition[startIndex]

	timeDiff <- stratColumn$modelTime - peakExtTimeMy
	index1 <- which(timeDiff == min(abs(timeDiff)))
	index2 <- which(timeDiff == -min(abs(timeDiff)))
	midIndex <- max(index1, index2)
	extMidHorizon <- stratColumn$stratPosition[midIndex]
		
	horizons <- list(extStart=extStartHorizon, extMid=extMidHorizon, extEnd=extEndHorizon)
	return(horizons)
}

# Test if mass extinction occurs
# Used by occurrencePlot() and fadLadPlot() to determine if an extinction box should be added
massExtinctionOccurred <- function(peakExtTimeMy, extDurationMy) {
	if ( !is.na(peakExtTimeMy) & !is.na(extDurationMy) ) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

# Add a box to a plot showing the time of extinction
# Used by occurrencePlot() and fadLadPlot()
addExtinctionBox <- function(stratColumn, peakExtTimeMy, extDurationMy, leftEdge, rightEdge, extinctionColor='red') {
	rgbColor <- col2rgb(extinctionColor)
	fillColor <- rgb(rgbColor[1], rgbColor[2], rgbColor[3], alpha=30, maxColorValue=255)	
	horizons <- extinctionHorizons(stratColumn=stratColumn, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)
	rect(leftEdge, horizons$extStart, rightEdge, horizons$extEnd, col=fillColor, border=fillColor)
	segments(leftEdge, horizons$extMid, rightEdge, horizons$extMid, col=extinctionColor, lwd=2)
}



