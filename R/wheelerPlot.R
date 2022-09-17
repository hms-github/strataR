#' @title Wheeler Plot
#'
#' @description Plot of regions of sedimentation, non-deposition, and erosion throug time.
#'
#' @details Plot is a chronostratigraphic diagram, also called a Wheeler plot, showing times and locations non-marine deposition (green), marine deposition (tan), subaerial erosion (the hiatus, black), and sediment starvation (white). Also shown are sediments that have been erosionally removed (the degradational vacuity, gray), as well as the location of the shore (blue). The combined hiatus and degradational vacuity represent the lacuna of an unconformity.
#'
#' @param basin an object of class [`basin`].
#' @param setting a string, either "valley" or "interfluve".
#'
#' @export
#'
#' @examples
#' data(sedBasin)
#' wheelerPlot(sedBasin, setting="valley")
#'

wheelerPlot <- function(basin, setting=c('valley', 'interfluve')) {
	setting <- match.arg(setting)
	marineColor <- 'tan'
	coastalPlainColor <- 'olivedrab3'
	minSediment <- 0.00001

	# Sediment history
	sedimentAccumulated <- basin$sedimentAccumulatedValley
	if (setting == 'interfluve') {
		sedimentAccumulated <- basin$sedimentAccumulatedInterfluve
	}

	# Start the plot
	plot(1, 1, type='n', xlim=c(basin$parameters$fallLineX, basin$parameters$marginWidth), ylim=c(0, basin$parameters$duration), las=1, xlab='Distance (km)', ylab='model time (m.y.)')

	# Marine sediments 
	# These are shown everywhere, will be masked by nonmarine, unconformity, and starvation
	graphics::rect(basin$parameters$fallLineX, 0, basin$parameters$marginWidth, basin$parameters$duration, col=marineColor, border=marineColor, lwd=0)

	# Overlay coastal plain sediments
	coastalPlainX <- c(basin$shore, 0, 0, basin$shore[1])
	last <- length(basin$timePoints)
	coastalPlainY <- c(basin$timePoints, basin$timePoints[last], basin$timePoints[1], basin$timePoints[1]) 
	graphics::polygon(coastalPlainX, coastalPlainY, col=coastalPlainColor, border=coastalPlainColor)
	
	# Add the subaerial degradational vacuity
	if (setting == 'valley') {
		deflection <- finalDeflection(basin)
		numPositions <- length(basin$positions)
		numTimes <- length(basin$timePoints)
		eroded <- matrix(data=FALSE, nrow=numTimes, ncol=numPositions)
		elevation <- basin$elevationProfileValley
		left <- basin$hiatusValley[ , 1:(numPositions-1)]
		right <- basin$hiatusValley[ , 2:numPositions]
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
	vacuityTimes <- basin$timePoints[vacuityCells[, 1]]
	vacuityPositions <- basin$positions[vacuityCells[, 2]]
	graphics::points(vacuityPositions, vacuityTimes, pch=16, cex=0.2, col='gray30')

	# Identify marine and nonmarine areas
	marine <- matrix(data=FALSE, nrow=nrow(sedimentAccumulated), ncol=ncol(sedimentAccumulated))
	rightEdge <- length(basin$positions)
	for (i in 1:length(basin$timePoints)) {
		shoreCell <- basin$shore[i] / basin$parameters$deltaX
		marine[i, shoreCell:rightEdge] <- TRUE
	}
	nonmarine <- !marine

	# Identify hiatus areas
	hiatus <- matrix(data=FALSE, nrow=nrow(sedimentAccumulated), ncol=ncol(sedimentAccumulated))
	hiatus[sedimentAccumulated < minSediment] <- TRUE  # any zero is a hiatus

	# Identify subaerial unconformity areas and plot them
	unconformity <- hiatus & nonmarine
	unconformityCells <- which(unconformity, arr.ind=TRUE)
	unconformityTimes <- basin$timePoints[unconformityCells[, 1]]
	unconformityPositions <- basin$positions[unconformityCells[, 2]]
	graphics::points(unconformityPositions, unconformityTimes, pch=16, cex=0.2, col='black')

	# Identify starvation areas and plot them
	starvation <- hiatus & marine
	starvationCells <- which(starvation, arr.ind=TRUE)
	starvationTimes <- basin$timePoints[starvationCells[, 1]]
	starvationPositions <- basin$positions[starvationCells[, 2]]
	graphics::points(starvationPositions, starvationTimes, pch=16, cex=0.2, col='white')
	
	# Add the shore
	graphics::points(basin$shore, basin$timePoints, type='l', lwd=2, col='dodgerblue')
}
