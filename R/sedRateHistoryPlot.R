#' @title Sedimentation Rate History Plot
#'
#' @description Plot the history of sedimentation rate in a stratigraphic column.
#'
#' @details Plots the sedimentation rates within a stratigraphic column, with sedimentation rates calculated over intervals of a given thickness.
#'
#' @param column an object of class [`column`].
#' @param stratBin a number (in m) indicating the stratigraphic thickness over which  
#'   average sedimentation rate should be calculated.
#' @param yAxisLabels a BOOL indicating whether labels should be included on the y-axis.
#' @param ... Additional graphical parameters to be passed to the plot.
#'
#' @export
#'
#' @examples
#' data(coluValley)
#' sedRateHistoryPlot(column=coluValley, stratBin=5)
#'

sedRateHistoryPlot <- function(column, stratBin=10, yAxisLabels=TRUE, ...) {
	rates <- sedRateAveraged(column, stratBin=stratBin)
	ylab <- c('Stratigraphic Position (m)')
	if (yAxisLabels == FALSE) {
		ylab = c('')
	}
	yLimits <- c(min(column$stratPosition), max(column$stratPosition))
	plot(rates$sedRate / 1000, rates$stratPosition, las=1, xlab='Sed. Rate (m/kyr)', ylab=ylab, ylim=yLimits, type='o', col='brown', pch=16, axes=yAxisLabels, frame.plot=yAxisLabels, ...)
	if (yAxisLabels == FALSE) {
		graphics::axis(1)
	}
}