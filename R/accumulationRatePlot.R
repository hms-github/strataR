#' @title Sediment Accumulation Rate Plot
#'
#' @description Cross-section of sedimentary basin showing the relative rates of sediment accumulation in nonmarine portions of the basin. 
#'
#' @details Relative rates of sediment accumulation are shown in gray scale, from white (lowest rates of accumulation) to black (highest rates). Rates are not indicated for marine portions of the basin, which is shown in uniform gray.
#'
#' @param basin an object of class [`basin`].
#' @param setting a string, either "valley" or "interfluve".
#' @param xlim,ylim Optional setting to set the range of the x-axis and y-axis that is 
#'   displayed.
#' @param ... Additional arguments to be passed to plot.
#'
#' @export
#'
#' @examples
#' data(sedBasin)
#' accumulationRatePlot(sedBasin, setting="valley")
#'

accumulationRatePlot <- function(basin, setting=c('valley', 'interfluve'), xlim=NULL, ylim=NULL, ...) {
	setting <- match.arg(setting)
		
	elevation <- 1
	sediment <- 1
	deflection <- finalDeflection(basin)
	if (setting == 'valley') {
	 	elevation <- basin$elevationProfileValley
	 	sediment <- basin$sedimentAccumulatedValley
	} else if (setting == 'interfluve') {
		elevation <- basin$elevationProfileInterfluve
	 	sediment <- basin$sedimentAccumulatedInterfluve
	} else {
	     warning("Setting must be 'valley' or 'interfluve'", call.=FALSE, immediate.=TRUE)
	}
	
	# Part 1: plot background of all sediments
	endIndex <- length(basin$timePoints)
	endTime <- basin$timePoints[endIndex]
	
	marineColor <- grDevices::gray(0.7)
	startSurface <- elevation[1, ] - deflection[1, ]
	endSurface <- elevation[endIndex, ]
	
	# handle any erosion of the starting surface at a valley
	if (setting == 'valley') {
		for (timeIndex in 2:endIndex) {
			overlyingSurface <- elevation[timeIndex, ]
			deflectedOverlyingSurface <- overlyingSurface - deflection[timeIndex, ]
			startSurface <- pmin(startSurface, deflectedOverlyingSurface)
		}
	}

	# If no specific range is set for x and y, create the default plot limits
	if (is.null(xlim)) {
		xlim=range(basin$positions)
	}
	if (is.null(ylim)) {
		ylim <- c(min(endSurface, startSurface), max(endSurface, startSurface))
	}
	
	plot(basin$positions, startSurface, type='n', xlim=xlim, ylim=ylim, xlab='Distance (km)', ylab='Elevation (m)', las=1)

	sedimentsY <- c(startSurface, rev(endSurface))
	sedimentsX <- c(basin$positions, rev(basin$positions))
	graphics::polygon(sedimentsX, sedimentsY, col=marineColor, lwd=0.5)
	
	# Part 2: plot non-marine accumulation rates
	finalElevations <- matrix(0, nrow=nrow(elevation), ncol=ncol(elevation))
	sediment[sediment < 0] <- 0
	marineAccumulation <- 9999   # will be replaced below
	numTimePoints <- length(basin$timePoints)
	
	# final elevation of base of strata
	finalElevations[1, ] <- elevation[1, ] - deflection[1, ]
	
	for (timeIndex in 2:numTimePoints) {                    # no accumulation at time 1
		surface <- elevation[timeIndex, ] - deflection[timeIndex, ]
		
		# save finalElevations to matrix
		finalElevations[timeIndex, ] <- surface
		
		# set marine values a constant
		marine <- which(elevation[timeIndex, ] < 0)
		sediment[timeIndex, marine] <- marineAccumulation
		
		if (setting == 'valley' & timeIndex<numTimePoints) {  
			# erode any underlying units (will just remove any unit even partly eroded)
			overlyingIndices <- seq(timeIndex+1, length(basin$timePoints))
			for (j in overlyingIndices) {
				deflectedOverlyingSurface <- elevation[j, ] - deflection[j, ]
				erosionPositions <- surface > deflectedOverlyingSurface
				sediment[timeIndex, erosionPositions] <- 0
			}
		}
	}
	
	# build x-y-z arrays for plotting
	vx <- c(matrix(basin$positions, nrow=nrow(sediment), ncol=ncol(sediment), byrow=TRUE))
	vy <- c(finalElevations)
	vz <- c(sediment)
	zero <- which(sediment==0)

	vx <- vx[-c(zero)]
	vy <- vy[-c(zero)]
	vz <- vz[-c(zero)]
	
	# rescaling because vz has a very odd distribution, strongly right-tailed without
	#   the marine values, plus a single spike for all the marine values
	nonmarinePoints <- (vz != marineAccumulation)
	vxNonmarine <- vx[nonmarinePoints]
	vyNonmarine <- vy[nonmarinePoints]
	vzNonmarine <- vz[nonmarinePoints]
	vzNonmarineRescaled <- 1 - rank(vzNonmarine)/length(vzNonmarine)
	
	
	graphics::points(vxNonmarine, vyNonmarine, col=grDevices::gray(vzNonmarineRescaled), pch=16, cex=0.5, lwd=3)
}
