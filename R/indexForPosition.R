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

indexForPosition <- function(positionKm, marginWidth, deltaX) {
	# positionKm, marginWidth, and deltaX all in km
	allPositions <- seq(0, marginWidth, deltaX)
	deviations <- abs(allPositions - positionKm)
	index <- which(deviations == min(deviations))
	index <- index[1]   # in case of a tie, choose the left-most position index
	return(index)
}