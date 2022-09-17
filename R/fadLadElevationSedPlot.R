#' @title Plot FADs or LADs, Elevation, and Sedimentation Rate
#'
#' @description Plot the number of FADs or LADs through a stratigraphic column, next to plots that show the change in elevation (meters above sea level) recorded in the column, and sedimentation rate.
#'
#' @details Plot consists of three panels, combining `fadLadPlot`, `elevationHistoryPlot`, and `sedRateHistoryPlot`. The left panel shows the number of first and last occurrences through a stratigraphic column. If a mass extinction was simulated, its peak and its duration are plotted as well. The sampling interval for binning first and last occurrences in specified. The middle panel shows the changes in elevation preserved through the stratigraphic column, in meters above sea level. The right panel shows sedimentation rate (in m / m.y.) through the stratigraphic column. The interval (in meters) over which sedimentation rates are calculated is specified.
#'
#' @param occurrences an object of class [`occurrences`].
#' @param column an object of class [`column`].
#' @param type a string of either "fad" or "lad" specifying which is to be plotted.
#' @param sampleSpacing a numeric value indicating the binning interval (in meters) for  
#'   tallying first or last occurrences in the stratigraphic column.
#' @param stratBin a numeric value indicating the stratigraphic interval (in meters) to be  
#'   binned for calculating sedimentation rates.
#' @param peakExtTimeMy a numeric value indicating when a mass extinction reached its peak  
#'   (in m.y.), if one was simulated. Use NA (the default) if no extinction was simulated.
#' @param extDurationMy a numeric value indicating when the duration of a mass extinction  
#'   (in m.y.), if one was simulated. Use NA (the default) if no extinction was simulated.
#'
#' @export
#'
#' @examples
#' data(occu)
#' data(coluValley)
#' fadLadElevationSedPlot(occurrences=occu, column=coluValley, type='lad', 
#'   sampleSpacing=0.5, stratBin=10)

#'

fadLadElevationSedPlot <- function(occurrences, column, type=c('fad', 'lad'), sampleSpacing=0.5, stratBin=10, peakExtTimeMy=NA, extDurationMy=NA) {
	opar <- graphics::par(no.readonly = TRUE)
	graphics::par(cex.axis=0.8)
	leftDivider <- 0.395  # picked by trial and error to give correct proportions
	rightDivider <- 0.675
	
	# left panel: first or last occurrences
	graphics::par(fig=c(0, leftDivider, 0, 1), mar=c(5, 4, 4, 0.5) + 0.1)
	fadLadPlot(type=type, occurrences=occurrences, column=column, sampleSpacing=sampleSpacing, xMax=NA, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)
	
	# middle panel: elevation
	graphics::par(fig=c(leftDivider, rightDivider, 0, 1), mar=c(5, 0, 4, 0.5) + 0.1, new=TRUE)
	elevationHistoryPlot(column, yAxisLabels=FALSE)
	
	# right panel: sedimentation rate
	graphics::par(fig=c(rightDivider, 1, 0, 1), mar=c(5, 0, 4, 2) + 0.1, new=TRUE)
	sedRateHistoryPlot(column, stratBin=stratBin, yAxisLabels= FALSE)
	
	graphics::par(opar)
}