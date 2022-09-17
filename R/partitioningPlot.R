#' @title Sediment Partitioning Plot
#' 
#' @description Plot showing how sediment is partitioned into nonmarine and marine areas through time.
#'
#' @details Plots shows the volume of sediment deposited in nonmarine vs. marine areas through time, as a fraction of the total amount of sediment. As the shore transgresses and regresses, owing to subsidence, eustasy, and sediment supply, deposition may occur preferentially in nonmarine or marine areas.
#'
#' @param basin an object of class [`basin`].
#'
#' @export
#'
#' @examples
#' data(sedBasin)
#' partitioningPlot(sedBasin)
#'

partitioningPlot <- function(basin) {
	marineColor <- 'tan'
	coastalPlainColor <- 'olivedrab3'
	timePoints <- basin$timePoints
	partition <- sedimentPartitioning(basin)
	
	plot(timePoints[-1], rep(1, length(timePoints[-1])), type='n', las=1, xlab='Model time (m.y.)', ylab='Proportion of sediment volume', ylim=c(0, 1))
	
	# Fill entire plot with marine color
	polygonX <- c(timePoints[-1], rev(timePoints[-1]))
	marinePolygonY <- c(rep(1, length(partition$totalVolume)), rep(0, length(partition$totalVolume)))
	graphics::polygon(polygonX, marinePolygonY, border=marineColor, col=marineColor)
	
	# Add nonmarine color on top
	correctedNonmarine <- partition$fractionNonmarine
	correctedNonmarine[is.na(correctedNonmarine)] <- 1.0
	nonmarinePolygonY <- c(correctedNonmarine, rep(0, length(partition$totalVolume)))
	graphics::polygon(polygonX, nonmarinePolygonY, border=coastalPlainColor, col=coastalPlainColor)
	
	# Add the partitioning line
	graphics::points(timePoints[-1], correctedNonmarine, type='l', lwd=1.5)
}
