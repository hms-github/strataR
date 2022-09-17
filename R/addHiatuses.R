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

addHiatuses <- function(basin, setting=c('valley', 'interfluve'), col='red', lwd=2, ...) {
	setting <- match.arg(setting)
	deflection <- finalDeflection(basin)
	numPositions <- length(basin$positions)
	
	hiatusX0 <- vector(mode='numeric')
	hiatusX1 <- vector(mode='numeric')
	hiatusY0 <- vector(mode='numeric')
	hiatusY1 <- vector(mode='numeric')
	
	if (setting == 'valley') {
	 	elevation <- basin$elevationProfileValley
	 	left <- basin$hiatusValley[ , 1:(numPositions-1)]
	 	right <- basin$hiatusValley[ , 2:numPositions]
	 	hiatusSegment <- left & right
	 	for (timeIndex in 2:(length(basin$timePoints)-1)) {
	 		surface <- elevation[timeIndex, ] - deflection[timeIndex, ]
			lowest <- timeIndex+1
			highest <- length(basin$timePoints)
			deflectedOverlyingSurfaces <- elevation[lowest:highest, ] - deflection[lowest:highest, ]
			lowestOverlyingSurface <- deflectedOverlyingSurfaces
			if (is.vector(lowestOverlyingSurface) != TRUE) {
				lowestOverlyingSurface <- apply(deflectedOverlyingSurfaces, MARGIN=2, FUN=min)
			}
			surface <- pmin(surface, lowestOverlyingSurface)
			origin <- which(hiatusSegment[timeIndex, ])
			final <- origin + 1
			hiatusX0 <- c(hiatusX0, basin$positions[origin])
			hiatusX1 <- c(hiatusX1, basin$positions[final])
			hiatusY0 <- c(hiatusY0, surface[origin])
			hiatusY1 <- c(hiatusY1, surface[final])
		}	
		
	 	segments(hiatusX0, hiatusY0, hiatusX1, hiatusY1, col=col, lwd=lwd, ...)
	} else if (setting == 'interfluve') {
		elevation <- basin$elevationProfileInterfluve
	 	left <- basin$hiatusInterfluve[ , 1:(numPositions-1)]
	 	right <- basin$hiatusInterfluve[ , 2:numPositions]
	 	hiatusSegment <- left & right
	 	for (timeIndex in 2:(length(basin$timePoints)-1)) {
	 		surface <- elevation[timeIndex, ] - deflection[timeIndex, ]
			lowest <- 1
			highest <- timeIndex-1
			deflectedUnderlyingSurfaces <- elevation[lowest:highest, ] - deflection[lowest:highest, ]
			highestUnderlyingSurface <- deflectedUnderlyingSurfaces
			if (is.vector(highestUnderlyingSurface) != TRUE) {
				highestUnderlyingSurface <- apply(deflectedUnderlyingSurfaces, MARGIN=2, FUN=max)
			}
			surface <- pmax(surface, highestUnderlyingSurface)
			origin <- which(hiatusSegment[timeIndex, ])
			final <- origin + 1
			hiatusX0 <- c(hiatusX0, basin$positions[origin])
			hiatusX1 <- c(hiatusX1, basin$positions[final])
			hiatusY0 <- c(hiatusY0, surface[origin])
			hiatusY1 <- c(hiatusY1, surface[final])
	 	}
	 	segments(hiatusX0, hiatusY0, hiatusX1, hiatusY1, col=col, lwd=lwd, ...)
	} else {
	     warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
	}
}