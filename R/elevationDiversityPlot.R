#' @title Elevation-Diversity Plot
#' 
#' @description Plot of diversity vs. elevation for species.
#'
#' @details Diagnostic plot to determine if diversity covaries with elevation; it shouldn't, but this might happen for stochastic reasons, especially if the number of species is small. Plot spans a specified range of elevations and records a species as present at a particular elevation if its probability of collection exceeds a minimum probability of collection. See `elevationDiversityResamplePlot()` for a version that adds confidence intervals.
#'
#' @param species an object of class [`species`].
#' @param minElevation,maxElevation a numeric value for the minimum and maximum elevation  
#'   to be considered, in meters above sea level.
#' @param pCrit a numeric value between 0 and 1 for the threshold probability of  
#'   collection for a species to be considered as occurring at a particular elevation.
#'
#' @export
#'
#' @examples
#' data(spec)
#' elevationDiversityPlot(species=spec, minElevation=20, maxElevation=200, pCrit=0.2)
#'

elevationDiversityPlot <- function(species, minElevation=0, maxElevation=150, pCrit=0.1) {
	numSpecies <- length(species$id)
	elevation <- seq(minElevation, maxElevation, 1)
	diversity <- rep(0, length(elevation))
	
	for (elev in 1:length(elevation)) {
		elevDiversity <- 0
		for (i in 1:numSpecies) {
			probCollection <- probCollection(PE=species$PE[i], ET=species$ET[i], PA=species$PA[i], elevation=elevation[elev])
			if (probCollection >= pCrit) {	                         # found
				elevDiversity <- elevDiversity + 1
			}
		}
		diversity[elev] <- elevDiversity
	}
	
	plot(elevation, diversity, ylim=c(0, max(diversity)), xlab='Elevation (m)', ylab='Richness (number of species)', type='l', lwd=2, las=1)
}