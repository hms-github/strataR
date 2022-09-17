# Occurrences Module

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
addExtinctionBox <- function(column, peakExtTimeMy, extDurationMy, leftEdge, rightEdge, extinctionColor='red') {
	rgbColor <- grDevices::col2rgb(extinctionColor)
	fillColor <- grDevices::rgb(rgbColor[1], rgbColor[2], rgbColor[3], alpha=30, maxColorValue=255)	
	horizons <- extinctionHorizons(stratColumn=column, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)
	graphics::rect(leftEdge, horizons$extStart, rightEdge, horizons$extEnd, col=fillColor, border=fillColor)
	graphics::segments(leftEdge, horizons$extMid, rightEdge, horizons$extMid, col=extinctionColor, lwd=2)
}


