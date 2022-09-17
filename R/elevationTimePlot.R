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

elevationTimePlot <- function(basin, locationKm, setting=c('valley', 'interfluve'), ...) {
	setting <- match.arg(setting)
	elevation <- elevationAtLocation(basin=basin, locationKm=locationKm, setting=setting)
	modelTime <- basin$timePoints
	plot(modelTime, elevation, type='l', lwd=2, col='brown', las=1, xlab='Model time (m.y.)', ylab='Elevation (m)', ...)
	abline(h=0, lty='dotted', col='gray')
}
