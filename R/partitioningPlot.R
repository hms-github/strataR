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

partitioningPlot <- function(basin) {
	marineColor <- 'tan'
	coastalPlainColor <- 'olivedrab3'
	timePoints <- basin$timePoints
	partition <- sedimentPartitioning(basin)
	
	plot(timePoints[-1], rep(1, length(timePoints[-1])), type='n', las=1, xlab='Model time (m.y.)', ylab='Proportion of sediment volume', ylim=c(0, 1))
	
	# Fill entire plot with marine color
	polygonX <- c(timePoints[-1], rev(timePoints[-1]))
	marinePolygonY <- c(rep(1, length(partition$totalVolume)), rep(0, length(partition$totalVolume)))
	polygon(polygonX, marinePolygonY, border=marineColor, col=marineColor)
	
	# Add nonmarine color on top
	correctedNonmarine <- partition$fractionNonmarine
	correctedNonmarine[is.na(correctedNonmarine)] <- 1.0
	nonmarinePolygonY <- c(correctedNonmarine, rep(0, length(partition$totalVolume)))
	polygon(polygonX, nonmarinePolygonY, border=coastalPlainColor, col=coastalPlainColor)
	
	# Add the partitioning line
	points(timePoints[-1], correctedNonmarine, type='l', lwd=1.5)
}
