#' @title Elevation History Plot
#' 
#' @description A plot of the history of elevation in a stratigraphic column.
#'
#' @details Plot shows how elevation above sea level (in m) preserved in a stratigraphic column changes within the column.
#'
#' @param column an object of class [`column`].
#' @param yAxisLabels a BOOL indicating whether labels should be included on the y-axis.
#' @param ... Additional graphical parameters to be passed to the plot.
#'
#' @export
#'
#' @examples
#' data(coluValley)
#' elevationHistoryPlot(coluValley)
#'

elevationHistoryPlot <- function(column, yAxisLabels=TRUE, ...) {
	ylab <- c('Stratigraphic Position (m)')
	if (yAxisLabels == FALSE) {
		ylab = c('')
	}
	yLimits <- c(min(column$stratPosition), max(column$stratPosition))
	plot(column$elevation, column$stratPosition, type='n', xlab='Elevation (m)', ylab=ylab, ylim=yLimits, las=1, axes=yAxisLabels, frame.plot=yAxisLabels, ...)
	if (yAxisLabels == FALSE) {
		graphics::axis(1)
	}
	graphics::abline(v=0, col='gray', lty='dotted')
	graphics::points(column$elevation, column$stratPosition, type='l', lwd=1.5)
}