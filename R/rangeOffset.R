#' @title Range Offset
#' 
#' @description Calculates range offset, the difference (in m.y.) between the time of origination and the age of the first occurrence, or between the time of extinction and the age of the last occurrence.
#'
#' @details The first and last occurrence of a fossil species seldom corresponds to the times of origination and extinction. The mismatch is called range offset.
#'
#' @param occurrences an object of class [`occurrences`].
#' @param species an object of class [`species`].
#'
#' @export
#' 
#' @return a data frame listing the species that occurred in a stratigraphic column, the stratigraphic positions (in meters) of their first and last occurrences, and the range offsets (in m.y.) of their first and last occurrences.
#'
#' @examples
#' data(occu)
#' data(spec)
#' rangeOffset(occurrences=occu, species=spec)
#'

rangeOffset <- function(occurrences, species) {
	# easier if occurrences and species are data frames, extracting only what is needed
	occurrencesDf <- data.frame(modelTime=occurrences$modelTime, stratPosition=occurrences$stratPosition, speciesId=occurrences$speciesId)
	speciesDf <- data.frame(id=species$id, origination=species$origination, extinction=species$extinction)
	
	occurringSpecies <- sort(unique(occurrencesDf$speciesId))
	numOccurringSpecies <- length(occurringSpecies)
	fad <- vector(mode='numeric', length=numOccurringSpecies)
	lad <- vector(mode='numeric', length=numOccurringSpecies)
	fadOffset <- vector(mode='numeric', length=numOccurringSpecies)
	ladOffset <- vector(mode='numeric', length=numOccurringSpecies)
	
	for (i in 1:numOccurringSpecies) {
		currentSpeciesId <- occurringSpecies[i]
		speciesOccurrence <- occurrencesDf[occurrencesDf$speciesId == currentSpeciesId, ]
		fad[i] <- min(speciesOccurrence$stratPosition)
		lad[i] <- max(speciesOccurrence$stratPosition)
		
		fadTime <- min(speciesOccurrence$modelTime)
		origination <- speciesDf$origination[currentSpeciesId]
		fadOffset[i] <- fadTime - origination
		
		ladTime <- max(speciesOccurrence$modelTime)
		extinction <- speciesDf$extinction[currentSpeciesId]
		ladOffset[i] <- extinction - ladTime
	}
	
	results <- data.frame(occurringSpecies, fad, lad, fadOffset, ladOffset)
	return(results)
}