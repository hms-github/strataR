#' @title Resampled Water Depth-Diversity Plot
#' 
#' @description A plot of diversity vs. water depth based on resampling.
#'
#' @details Diagnostic plot to determine if diversity covaries with water depth; it shouldn't, but this might happen for stochastic reasons, especially if the number of species is small. Plot spans a specified range of water depths and resamples to test if each species is present at a particular water depth. See `marineWaterDepthDiversityPlot()` for a faster version that lacks confidence intervals and is based on a fixed threshold probability of collection.
#'
#' @param marineSpecies an object of class [`marineSpecies`].
#' @param minWaterDepth,maxWaterDepth a numeric value for the minimum and maximum water  
#'   depth to be considered, in meters below sea level.
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
#' data(marspec)
#' marineWaterDepthDiversityResamplePlot(marineSpecies=marspec, minWaterDepth=0, 
#'   maxWaterDepth=150, numTrials=100)
#'

marineWaterDepthDiversityResamplePlot <- function(marineSpecies, minWaterDepth=0, maxWaterDepth=150, time=0, numTrials=100, conf=0.50) {
	# easier if species is a data frame
	speciesDf <- data.frame(id=marineSpecies$id, ancestor=marineSpecies$ancestor, origination=marineSpecies$origination, extinction=marineSpecies$extinction, PD=marineSpecies$PD, DT=marineSpecies$DT, PA=marineSpecies$PA)
	
	selectedSpecies <- speciesDf[which(speciesDf$origination <= time & speciesDf$extinction >= time), ]
	numSpecies <- nrow(selectedSpecies)
	
	waterDepth <- seq(minWaterDepth, maxWaterDepth, 1)
	diversityMid <- rep(0, length(waterDepth)) # median
	diversityLow <- rep(0, length(waterDepth)) # lower bound
	diversityHigh <- rep(0, length(waterDepth)) # higher bound
	significance <- 1 - conf
	lowBound <- significance/2
	highBound <- 1 - lowBound
	
	for (depth in 1:length(waterDepth)) {
		trialDiversity <- rep(0, numTrials)
		for (trial in 1:numTrials) {
			iterationDiversity <- 0
			for (i in 1:numSpecies) {
				random <- stats::runif(1)
				marineProbCollection <- marineProbCollection(PD=selectedSpecies$PD[i], DT=selectedSpecies$DT[i], PA=selectedSpecies$PA[i], waterDepth=waterDepth[depth])
				if (random <= marineProbCollection) {	                         # found
					iterationDiversity <- iterationDiversity + 1
				}
			}
			trialDiversity[trial] <- iterationDiversity
		}
		diversityMid[depth] <- stats::quantile(trialDiversity, 0.50)
		diversityLow[depth] <- stats::quantile(trialDiversity, lowBound)
		diversityHigh[depth] <- stats::quantile(trialDiversity, highBound)
	}

	yLimits <- c(0, max(diversityHigh))
	plot(waterDepth, diversityMid, ylim=yLimits, xlab='Water Depth (m)', ylab='Richness (number of species)', type='n', lwd=2, las=1)

	polyX <- c(waterDepth, rev(waterDepth))
	polyY <- c(diversityLow, rev(diversityHigh))
	graphics::polygon(polyX, polyY, col='lightgray', border='lightgray')
	
	graphics::points(waterDepth, diversityMid, type='l', lwd=2, col='black')
}