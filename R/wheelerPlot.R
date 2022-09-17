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

wheelerPlot <- function(x, setting=c('valley', 'interfluve')) {
	setting <- match.arg(setting)
	marineColor <- 'tan'
	coastalPlainColor <- 'olivedrab3'
	minSediment <- 0.00001

	# Sediment history
	sedimentAccumulated <- x$sedimentAccumulatedValley
	if (setting == 'interfluve') {
		sedimentAccumulated <- x$sedimentAccumulatedInterfluve
	}

	# Start the plot
	dev.new()
	plot(1, 1, type='n', xlim=c(x$parameters$fallLineX, x$parameters$marginWidth), ylim=c(0, x$parameters$duration), las=1, xlab='Distance (km)', ylab='model time (m.y.)')

	# Marine sediments 
	# These are shown everywhere, will be masked by nonmarine, unconformity, and starvation
	rect(x$parameters$fallLineX, 0, x$parameters$marginWidth, x$parameters$duration, col=marineColor, border=marineColor, lwd=0)

	# Overlay coastal plain sediments
	coastalPlainX <- c(x$shore, 0, 0, x$shore[1])
	last <- length(x$timePoints)
	coastalPlainY <- c(x$timePoints, x$timePoints[last], x$timePoints[1], x$timePoints[1]) 
	polygon(coastalPlainX, coastalPlainY, col=coastalPlainColor, border=coastalPlainColor)
	
	# Add the subaerial degradational vacuity
	if (setting == 'valley') {
		deflection <- finalDeflection(x)
		numPositions <- length(x$positions)
		numTimes <- length(x$timePoints)
		eroded <- matrix(data=FALSE, nrow=numTimes, ncol=numPositions)
		elevation <- x$elevationProfileValley
		left <- x$hiatusValley[ , 1:(numPositions-1)]
		right <- x$hiatusValley[ , 2:numPositions]
		hiatusSegment <- left & right
		for (timeIndex in 2:(numTimes-1)) {
			surface <- elevation[timeIndex, ] - deflection[timeIndex, ]
			lowest <- timeIndex+1
			highest <- numTimes
			deflectedOverlyingSurfaces <- elevation[lowest:highest, ] - deflection[lowest:highest, ]
			lowestOverlyingSurface <- deflectedOverlyingSurfaces
			if (is.vector(lowestOverlyingSurface) != TRUE) {
				lowestOverlyingSurface <- apply(deflectedOverlyingSurfaces, MARGIN=2, FUN=min)
			}
			eroded[timeIndex, ] <- (surface - lowestOverlyingSurface) > 1
		}
	}
	vacuityCells <- which(eroded, arr.ind=TRUE)
	vacuityTimes <- x$timePoints[vacuityCells[, 1]]
	vacuityPositions <- x$positions[vacuityCells[, 2]]
	points(vacuityPositions, vacuityTimes, pch=16, cex=0.2, col='gray30')

	# Identify marine and nonmarine areas
	marine <- matrix(data=FALSE, nrow=nrow(sedimentAccumulated), ncol=ncol(sedimentAccumulated))
	rightEdge <- length(x$positions)
	for (i in 1:length(x$timePoints)) {
		shoreCell <- x$shore[i] / x$parameters$deltaX
		marine[i, shoreCell:rightEdge] <- TRUE
	}
	nonmarine <- !marine

	# Identify hiatus areas
	hiatus <- matrix(data=FALSE, nrow=nrow(sedimentAccumulated), ncol=ncol(sedimentAccumulated))
	hiatus[sedimentAccumulated < minSediment] <- TRUE  # any zero is a hiatus

	# Identify subaerial unconformity areas and plot them
	unconformity <- hiatus & nonmarine
	unconformityCells <- which(unconformity, arr.ind=TRUE)
	unconformityTimes <- x$timePoints[unconformityCells[, 1]]
	unconformityPositions <- x$positions[unconformityCells[, 2]]
	points(unconformityPositions, unconformityTimes, pch=16, cex=0.2, col='black')

	# Identify starvation areas and plot them
	starvation <- hiatus & marine
	starvationCells <- which(starvation, arr.ind=TRUE)
	starvationTimes <- x$timePoints[starvationCells[, 1]]
	starvationPositions <- x$positions[starvationCells[, 2]]
	points(starvationPositions, starvationTimes, pch=16, cex=0.2, col='white')
	
	# Add the shore
	points(x$shore, x$timePoints, type='l', lwd=2, col='dodgerblue')
}
