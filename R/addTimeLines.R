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

addTimeLines <- function(basin, setting=c('valley', 'interfluve'), timeLines, ...) {
	setting <- match.arg(setting)
	elevation <- 1
	deflection <- finalDeflection(basin)
	
	if (setting == 'valley') {
	 	elevation <- basin$elevationProfileValley
	} else if (setting == 'interfluve') {
		elevation <- basin$elevationProfileInterfluve
	} else {
	     warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
	}
	
	for (i in 1:length(timeLines)) {
		timeIndex <- which(abs(basin$timePoints - timeLines[i]) == min(abs(basin$timePoints - timeLines[i])))
		surface <- elevation[timeIndex, ] - deflection[timeIndex, ]
		
		# in the valley where erosion can occur, correct all older positions to be no higher than the current profile. On the interfluve, ensure that the profile is no lower than any previous profile.
		if (setting == 'valley') {
			overlyingIndices <- seq(timeIndex+1, length(basin$timePoints))
			for (j in overlyingIndices) {
				deflectedOverlyingSurface <- elevation[j, ] - deflection[j, ]
				surface <- pmin(surface, deflectedOverlyingSurface)
			}
		} else if (setting == 'interfluve') {
			underlyingIndices <- seq(1, timeIndex-1)
			for (j in underlyingIndices) {
				deflectedUnderlyingSurface <- elevation[j, ] - deflection[j, ]
				surface <- pmax(surface, deflectedUnderlyingSurface)
			}
		} else {
			warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
		}
		
		points(basin$positions, surface, type='l', ...)
	}
}
