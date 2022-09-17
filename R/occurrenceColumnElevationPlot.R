#' @title Plot Fossil Occurrences and Elevation
#' 
#' @description Plot a stratigraphic column, a chart showing fossil occurrences and ranges, and elevation through the column.
#'
#' @details Produces a three panel plot. The left plot is the stratigraphic column, produced with `plot.column`. Yellow indicates channel deposits, green indicates floodplain deposits, and tan corresponds to marine deposits. Hiatuses are shown with a red line. 
#' The middle plot shows fossil occurrences and ranges, produced with `plot.occurrences`. Occurrences are shown with solid circles, except for singletons (species occurring only once), which are shown with open circles. Crosses indicate horizons that correspond to the times of origination and extinction. Blue lines indicate the preserved range of a species (i.e., connecting the first and last occurrence), and gray lines indicate the interval of which a species was extant (i.e., between the horizons corresponding to the times of origination and extinction. If a mass extinction was simulated, the stratigraphic interval corresponding to the time of extinction and its peak will be indicated on the stratigraphic column.
#' The right plot shows the changes in elevation (meters above sea level) preserved in the column. It is generated with `elevationHistoryPlot`.
#'
#' @param occurrences an object of class [`occurrences`].
#' @param column an object of class [`column`].
#' @param species an object of class [`species`].
#' @param peakExtTimeMy a numeric value indicating when a mass extinction reached its peak  
#'   (in m.y.), if one was simulated. Use NA (the default) if no extinction was simulated.
#' @param extDurationMy a numeric value indicating when the duration of a mass extinction  
#'   (in m.y.), if one was simulated. Use NA (the default) if no extinction was simulated.
#' @param orderedByLad a boolean indicating how species will be sorted on the occurrence  
#'   chart. If TRUE, they are ordered by their last occurrences; if FALSE, they are  
#'   ordered by their first occurrences.
#'
#' @export
#' 
#' @examples
#' data(occu)
#' data(coluValley)
#' data(spec)
#' occurrenceColumnElevationPlot(occurrences=occu, column=coluValley, species=spec)
#'

occurrenceColumnElevationPlot <- function(occurrences, column, species, peakExtTimeMy=NA, extDurationMy=NA, orderedByLad=TRUE) {
	opar <- graphics::par(no.readonly = TRUE)
	leftDivider <- 0.20   # separates ranges from elevation plot
	rightDivider <- 0.80  # separates elevation plot from strat column
	
	# Left panel: stratigraphic column
	graphics::par(fig=c(0, leftDivider, 0, 1), mar=c(5, 4, 4, 0) + 0.1)
	plot(column)
	
	# Middle panel: Range chart
	graphics::par(fig=c(leftDivider, rightDivider, 0, 1), mar=c(5, 0, 4, 0) + 0.1, new=TRUE)
	plot(occurrences, column=column, species=species, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy, orderedByLad=orderedByLad, yAxisLabels=FALSE)
	
	# Right panel: elevation plot
	graphics::par(fig=c(rightDivider, 1, 0, 1), mar=c(5, 0, 4, 2) + 0.1, new=TRUE)
	elevationHistoryPlot(column, yAxisLabels=FALSE)
	
	graphics::par(opar)
}
