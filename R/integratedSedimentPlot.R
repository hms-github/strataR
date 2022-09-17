#' @title Integrated Sediment Plot
#'
#' @description Plot of the total sediment deposited through time.
#'
#' @details Plot shows the actual volume of sediment deposited through time versus the targeted volume of sediment. 
#'
#' @param basin an object of class [`basin`].
#'
#' @export
#'
#' @examples
#' data(sedBasin)
#' integratedSedimentPlot(sedBasin)
#'

integratedSedimentPlot <- function(basin) {
	plot(basin$timePoints[-1], basin$integratedSediment[-1], type='o', pch=16, cex=0.5, las=1, xlab='Model time (m.y.)', ylab='Integrated sediment volume', ylim=c(0, max(basin$integratedSediment)))
	graphics::points(basin$timePoints, basin$targetSediment, type='l', col='brown', lty='dashed', lwd=2)
	labelX <- 0.30 * max(basin$timePoints)
	labelY <- 0.30 * max(basin$integratedSediment)
	graphics::text(labelX, labelY, 'target sediment volume', pos=3, col='brown')
	labelY <- 0.25 * max(basin$integratedSediment)
	graphics::text(labelX, labelY, 'actual sediment volume', pos=3, col='black')
}
