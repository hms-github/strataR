#' @title First and Last Occurrences (FADs and LADs)
#'
#' @description `fadLad` returns the first and last occurrence of every species in a column.
#'
#' @details Finds the first and last occurrence of every species in a stratigraphic column, reporting the stratigraphic position (in m) and the age of that horizon (in m.y.). Also reports the number of occurrences of each species.
#'
#' @param occurrences an object of class [`occurrences`].
#'
#' @export
#' 
#' @return A data frame including the species ID of those that occur in a column, the stratigraphic position (in m) of their FAD and LAD, the model times (in m.y.) of the FAD and LAD, and the number of occurrences in the column. 
#'
#' @examples
#' data(occu)
#' fadLad(occurrences=occu)
#'

fadLad <- function(occurrences) {
	# convert occurrences from class 'occurrences' to a data.frame
	# extracting only what is needed
	occdf <- data.frame(modelTime=occurrences$modelTime, stratPosition=occurrences$stratPosition, speciesId=occurrences$speciesId)
	
	occurringSpecies <- sort(unique(occdf$speciesId))
	fad <- vector(mode='numeric', length=length(occurringSpecies))
	lad <- vector(mode='numeric', length=length(occurringSpecies))
	fadTime <- vector(mode='numeric', length=length(occurringSpecies))
	ladTime <- vector(mode='numeric', length=length(occurringSpecies))
	numOccurrences <- vector(mode='numeric', length=length(occurringSpecies))
	
	for (i in 1:length(occurringSpecies)) {
		speciesOccurrence <- occdf[occdf$speciesId == occurringSpecies[i], ]
		fad[i] <- min(speciesOccurrence$stratPosition)
		lad[i] <- max(speciesOccurrence$stratPosition)
		fadTime[i] <- min(speciesOccurrence$modelTime)
		ladTime[i] <- max(speciesOccurrence$modelTime)
		numOccurrences[i] <- length(speciesOccurrence$stratPosition)
	}
	
	fadsLads <- as.data.frame(cbind(occurringSpecies, fad, lad, fadTime, ladTime, numOccurrences))
	return(fadsLads)
}