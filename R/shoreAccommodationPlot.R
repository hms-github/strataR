#' Create a ternary plot
#'
#' Creates a ternary plot from three variables that sum to 100%, such as petrographic data in geology. Specifies the labels for the three apexes of the triangle, whether a grid should be shown, and the spacing of the grid lines
#'
#' @param \code{TRUE}
#'
#' @export
#' 
#' @return \code{TRUE}
#'
#' @examples
#' ternaryPlot(myData, labels=("Q", "F", "L"))
#'

shoreAccommodationPlot <- function(basin, shoreRates) {
	shoreRates <- shoreAccommodationRates(basin)
	minRate <- min(0, shoreRates$accommodation)
	maxRate <- max(0, shoreRates$accommodation)
	plot(shoreRates$modelTime, shoreRates$accommodation, type='l', xlab='Model time (m.y.)', ylab='Accommodation rate (m/m.y.)', ylim=c(minRate, maxRate), las=1, lwd=2)
	
	# Zero-rate line
	abline(h=0, col='black')
	
	# Subsidence rates at the shore
	points(shoreRates$modelTime, shoreRates$shoreSubsidence, type='l', col='gray', lty='dotted')
	
	# Times in which the shore reverses its direction of movement
	# Note that this only finds times when the shore position is not changing
	# If it changes immediately from regression to transgression, or vice versa, without pausing for a time of no change, this will not find that reversal
	shoreMoveSign <- shoreRates$shoreMovement / abs(shoreRates$shoreMovement)
	shoreMoveSign[is.nan(shoreMoveSign)] <- 0
	abline(v=which(shoreMoveSign == 0)*basin$parameters$timeStep, col='gray', lty='dashed')
}