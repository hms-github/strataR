#' @title Water Depth-Diversity Plot
#' 
#' @description Plot of diversity vs. water depth for marine species.
#'
#' @details Diagnostic plot to determine if diversity covaries with water depth; it shouldn't, but this might happen for stochastic reasons, especially if the number of species is small. Plot spans a specified range of water depths and records a species as present at a particular water depth if its probability of collection exceeds a minimum probability of collection. See `marineWaterDepthDiversityResamplePlot()` for a version that adds confidence intervals.
#' 
#' @param marineSpecies an object of class [`marineSpecies`].
#' @param minWaterDepth,maxWaterDepth a numeric value for the minimum water depth to be  
#'   considered, in meters below sea level.
#' @param pCrit a numeric value between 0 and 1 for the threshold probability of  
#'   collection for a species to be considered as occurring at a particular water depth.
#'
#' @export
#'
#' @examples
#' data(marspec)
#' marineWaterDepthDiversityPlot(marineSpecies=marspec, minWaterDepth=0, 
#'   maxWaterDepth=150, pCrit=0.1)
#'

marineWaterDepthDiversityPlot <- function(marineSpecies, minWaterDepth=0, maxWaterDepth=150, pCrit=0.1) {
	numSpecies <- length(marineSpecies$id)
	waterDepth <- seq(minWaterDepth, maxWaterDepth, 1)
	diversity <- rep(0, length(waterDepth))
	
	for (depth in 1:length(waterDepth)) {
		waterDepthDiversity <- 0
		for (i in 1:numSpecies) {
			marineProbCollection <- marineProbCollection(PD=marineSpecies$PD[i], DT=marineSpecies$DT[i], PA=marineSpecies$PA[i], waterDepth=waterDepth[depth])
			if (marineProbCollection >= pCrit) {	                         # found
				waterDepthDiversity <- waterDepthDiversity + 1
			}
		}
		diversity[depth] <- waterDepthDiversity
	}
	
	plot(waterDepth, diversity, ylim=c(0, max(diversity)), xlab='Water Depth (m)', ylab='Richness (number of species)', type='l', lwd=2, las=1)
}