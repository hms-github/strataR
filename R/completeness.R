#' @title Completeness
#'
#' @description Completeness of the fossil record.
#'
#' @details `completeness` calculates the percentage of species that were preserved at least once in a column (completeness), the percentage of species that had only two occurrences in a column (two-timers), and the percentage of species that occurred only once (singletons).
#'
#' @param occurrences an object of class [`occurrences`].
#' @param species an object of class [`species`].
#'
#' @export
#' 
#' @return A list of the total number species, the number of species that occurred at least once, the number species that occurred only once, and the number of species that occurred only twice.
#'
#' @examples
#' data(occu)
#' data(spec)
#' completeness(occurrences=occu, species=spec)
#'

completeness <- function(occurrences, species) {
	totalSpecies <- length(species)	
	tallies <- table(occurrences$speciesId)
	occurringSpecies <- length(tallies)
	singletons <- length(tallies[tallies == 1])
	twoTimers  <- length(tallies[tallies == 2])
	
	percentOccurrence <- round(occurringSpecies/totalSpecies*100, 1)
	percentSingletons <- round(singletons/totalSpecies*100, 1)
	percentTwoTimers  <- round(twoTimers/totalSpecies*100, 1)
	
	cat("There are ", totalSpecies, " species in the simulation, of which\n", sep='')
	cat("  ", occurringSpecies, " (", percentOccurrence, "%) left a fossil record, \n", sep='')
	cat("  ", singletons, " (", percentSingletons, "%) occur only once, and\n", sep='')
	cat("  ", twoTimers, " (", percentTwoTimers, "%) occur only twice.\n", sep='')
	cat(" \n")
	
	results <- list(totalSpecies=totalSpecies, occurringSpecies=occurringSpecies, singletons=singletons, twoTimers=twoTimers)
	return(results)
}
