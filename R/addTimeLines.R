#' @title Time Lines
#'
#' @description Adds time lines to a basin plot.
#'
#' @details Time lines of any desired spacing are added to a basin plot.
#'
#' @param basin an object of class [`basin`].
#' @param setting a string, either "valley" or "interfluve".
#' @param timeLines a vector of times (in m.y.) for which time lines should be added.
#' @param ... Optional arguments to style how the time lines are displayed.
#'
#' @export
#' 
#' @examples
#' data(sedBasin)
#' plot(sedBasin)
#' addTimeLines(sedBasin, setting='valley', timeLines=seq(0.1, 2.9, 0.1))
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
		
		graphics::points(basin$positions, surface, type='l', ...)
	}
}
