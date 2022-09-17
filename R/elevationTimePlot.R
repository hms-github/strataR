#' @title Elevation–Time Plot
#'
#' @description Plot of elevation through time.
#'
#' @details Plot showing how the elevation (in m above sea level) at a location in the basin changes over time. 
#'
#' @param basin an object of class [`basin`].
#' @param locationKm the location in km as a distance from the left edge of the basin
#' @param setting a string, either "valley" or "interfluve".
#' @param ... Optional arguments to be passed to the plot.
#'
#' @export
#'
#' @examples
#' data(sedBasin)
#' elevationTimePlot(sedBasin, locationKm=100, setting="valley")
#'

elevationTimePlot <- function(basin, locationKm, setting=c('valley', 'interfluve'), ...) {
	setting <- match.arg(setting)
	elevation <- elevationAtLocation(basin=basin, locationKm=locationKm, setting=setting)
	modelTime <- basin$timePoints
	plot(modelTime, elevation, type='l', lwd=2, col='brown', las=1, xlab='Model time (m.y.)', ylab='Elevation (m)', ...)
	graphics::abline(h=0, lty='dotted', col='gray')
}
