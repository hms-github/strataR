#' @title Shoreline Accommodation Plot
#'
#' @description Plot of the rate of accommodation at the shoreline through time.
#'
#' @details Shows the rate of accommodation (subsidence plus sea-level rise) at the shoreline through time. This will vary as rates of accommodation and eustasy vary through time, as well as with lateral shifts in the shoreline.
#'
#' @param basin an object of class [`basin`].
#'
#' @export
#'
#' @examples
#' data(sedBasin)
#' shoreAccommodationPlot(sedBasin)
#'

shoreAccommodationPlot <- function(basin) {
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