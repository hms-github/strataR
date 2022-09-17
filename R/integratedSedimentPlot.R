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

integratedSedimentPlot <- function(basin) {
	plot(basin$timePoints[-1], basin$integratedSediment[-1], type='o', pch=16, cex=0.5, las=1, xlab='Model time (m.y.)', ylab='Integrated sediment volume', ylim=c(0, max(basin$integratedSediment)))
	points(basin$timePoints, basin$targetSediment, type='l', col='brown', lty='dashed', lwd=2)
	labelX <- 0.30 * max(basin$timePoints)
	labelY <- 0.30 * max(basin$integratedSediment)
	text(labelX, labelY, 'target sediment volume', pos=3, col='brown')
	labelY <- 0.25 * max(basin$integratedSediment)
	text(labelX, labelY, 'actual sediment volume', pos=3, col='black')
}
