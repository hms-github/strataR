#' @title Resampled Elevation-Diversity Plot
#' 
#' @description A plot of diversity vs. elevation based on resampling.
#'
#' @details Diagnostic plot to determine if diversity covaries with elevation; it shouldn't, but this might happen for stochastic reasons, especially if the number of species is small. Plot spans a specified range of elevations and resamples to test if each species is present at a particular elevation. See `elevationDiversityPlot()` for a faster version that lacks confidence intervals and is based on a fixed threshold probability of collection.
#'
#' @param species an object of class [`species`].
#' @param minElevation,maxElevation a numeric value for the minimum and maximum elevation  
#'   to be considered, in meters above sea level.
#' @param time a numeric value, in m.y., that limits analysis to species that originated  
#'   before a specified time and that went extinct after that same time. Generally left at  
#'   default to include species alive at beginning of simulation.
#' @param numTrials a numeric value for the number of resampling trials.
#' @param conf is a numeric value (0 to 1) indicating the confidence level used for  
#'   drawing the uncertainty envelope
#'
#' @export
#'
#' @examples
#' data(spec)
#' elevationDiversityResamplePlot(species=spec, minElevation=0, maxElevation=150, 
#'   numTrials=100) 
#'

elevationDiversityResamplePlot <- function(species, minElevation=0, maxElevation=150, time=0, numTrials=100, conf=0.50) {
	# easier if species is a data frame
	speciesDf <- data.frame(id=species$id, ancestor=species$ancestor, origination=species$origination, extinction=species$extinction, PE=species$PE, ET=species$ET, PA=species$PA, Aff=species$Aff)
	
	selectedSpecies <- speciesDf[which(speciesDf$origination <= time & speciesDf$extinction >= time), ]
	numSpecies <- nrow(selectedSpecies)
	
	elevation <- seq(minElevation, maxElevation, 1)
	diversityMid <- rep(0, length(elevation)) # median
	diversityLow <- rep(0, length(elevation)) # lower bound
	diversityHigh <- rep(0, length(elevation)) # higher bound
	significance <- 1 - conf
	lowBound <- significance/2
	highBound <- 1 - lowBound
	
	for (elev in 1:length(elevation)) {
		trialDiversity <- rep(0, numTrials)
		for (trial in 1:numTrials) {
			iterationDiversity <- 0
			for (i in 1:numSpecies) {
				random <- stats::runif(1)
				probCollection <- probCollection(PE=selectedSpecies$PE[i], ET=selectedSpecies$ET[i], PA=selectedSpecies$PA[i], elevation=elevation[elev])
				if (random <= probCollection) {	                         # found
					iterationDiversity <- iterationDiversity + 1
				}
			}
			trialDiversity[trial] <- iterationDiversity
		}
		diversityMid[elev] <- stats::quantile(trialDiversity, 0.50)
		diversityLow[elev] <- stats::quantile(trialDiversity, lowBound)
		diversityHigh[elev] <- stats::quantile(trialDiversity, highBound)
	}

	yLimits <- c(0, max(diversityHigh))
	plot(elevation, diversityMid, ylim=yLimits, xlab='Elevation (m)', ylab='Richness (number of species)', type='n', lwd=2, las=1)

	polyX <- c(elevation, rev(elevation))
	polyY <- c(diversityLow, rev(diversityHigh))
	graphics::polygon(polyX, polyY, col='lightgray', border='lightgray')
	
	graphics::points(elevation, diversityMid, type='l', lwd=2, col='black')
}