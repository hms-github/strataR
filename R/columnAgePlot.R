#' @title Column-age Plot
#' 
#' @description A plot of model age through a stratigraphic column (a time-depth plot).
#'
#' @details Plot shows the age of each horizon within a column.
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

columnAgePlot <- function(column, yAxisLabels=TRUE, ...) {
	ylab <- c('Stratigraphic Position (m)')
	if (yAxisLabels == FALSE) {
		ylab = c('')
	}
	yLimits <- c(min(column$stratPosition), max(column$stratPosition))
	plot(column$modelTime, column$stratPosition, type='n', xlab='Model time (m.y.)', ylab=ylab, ylim=yLimits, las=1, axes=yAxisLabels, frame.plot=yAxisLabels, ...)
	if (yAxisLabels == FALSE) {
		graphics::axis(1)
	}
	graphics::abline(v=0, col='gray', lty='dotted')
	graphics::points(column$modelTime, column$stratPosition, type='l', lwd=1.5)
}