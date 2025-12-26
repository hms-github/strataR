#' @title Stratigraphic Columns
#'
#' @description Create a stratigraphic column from a basin.
#'
#' @details Create a stratigraphic column from anywhere in the basin by specifying the location of the column (the distance from the left edge of the basin) and by specifying whether the column is located along a valley or an interfluve. The probability of a channel occurring in any time step and the depth of the channel control control the relative proportions of chanel and floodplain facies.
#'
#' @param basin an object of class [`basin`].
#' @param locationKm a number indicating where in the basin the column will be located,  
#'   specified as the distance (in km) from the left edge of the basin.
#' @param setting a string, either "valley" or "interfluve", specifying where the colulmn  
#'   will be located.
#' @param pChannel a number from 0 to 1 specifying the probability that a channel will be  
#'   present at any time step. Floodplain facies are deposited if a channel is not  
#'   present, so 1-pChannel can be considered the probability of floodplain facies.
#' @param channelDepth a number specifying the depth of a channel, in m.
#' @param x,object an object of class `column`.
#' @param stratRange a vector specifying (in m) the range of y-values on the column plot.
#' @param axes a boolean specifying whether to display the axes on the column plot.
#' @param ... additional arguments to be passed to the column plot.
#'
#' @export
#' 
#' @return `column` returns an object of class "column", which includes print and summary methods.
#' 
#' A column object consists of a data frame with five columns; each row corresponds to one sedimentary deposit. `modelTime` is the model time step in millions of years, `thickness` is the thickness of the deposit in meters, `stratPosition` is the height (in meters) within the stratigraphic column of the base of the deposit, `facies` is the facies of the deposit (floodplain, channel, marine, or carbonate), and `elevation` is the elevation (in meters above sea level) at which the deposit accumulated. Water depths therefore have negative values.
#'
#' @examples
#' data(sedBasin)
#' colu <- column(sedBasin, locationKm=200, setting='valley', pChannel=0.1, channelDepth=2)
#' summary(colu)
#' plot(colu)
#' 
#' @rdname column
#' @export column

column <- function(basin, locationKm, setting=c('valley', 'interfluve'), pChannel, channelDepth) {	
	setting <- match.arg(setting)
	if (setting != 'valley' && setting != 'interfluve') {
	     warning("Setting must be 'valley' or 'interfluve'")
	}
	
	elevation <- elevationAtLocation(basin, locationKm, setting)             # m
	sedimentation <- sedimentationAtLocation(basin, locationKm, setting)     # m
	timeStep <- basin$parameters$timeStep                                    # m.y.
	
	modelTime     <- basin$timePoints
	thickness     <- rep( 0, length(modelTime))
	stratPosition <- rep( 0, length(modelTime))
	facies        <- rep('', length(modelTime))
	column <- data.frame(sedimentation, thickness, stratPosition, facies)
		
	# Build column
	for (i in 2:length(modelTime)) { # skip the first time step, the base of column
		if (sedimentation[i] > 0) {					# accumulation
			if (elevation[i] < 0) {							# marine
				column <- depositMarine(column=column, timeStep=i)
			} else if (stats::runif(1) > pChannel) {				# floodplain
				column <- depositFloodplain(column, i)
			} else {										# channel
				if (channelDepth <= sedimentation[i]) {				# aggrading
					column <- depositNonincisingChannel(column=column, timeStep=i)
				} else {											# incising
					column <- depositIncisingChannel(column=column, timeStep=i, channelDepth=channelDepth)
				}
			}          
		} else if (sedimentation[i] == 0) {			# bypass 
			column <- nondepositionDuringBypass(column=column, timeStep=i)
		} else {									# erosion
			column <- erodeAtUnconformity(column=column, timeStep=i)
		}
		
		# Test that column matches the basin
		tolerance <- 0.001
		if (column$stratPosition[i] > sum(column$sedimentation[1:i]) + tolerance) {
			message <- paste("Warning column is too thick in time ", i, ", column: ", column$stratPosition[i], ", basin: ", sum(column$sedimentation[1:i]))
			print(message)
		}
		if (column$stratPosition[i] < sum(column$sedimentation[1:i]) - tolerance) {
			message <- paste("Warning column is too thin in time ", i, ", column: ", column$stratPosition[i], ", basin: ",  sum(column$sedimentation[1:i]))
			print(message)
		}
	}

	# If necessary, explain why column might be slightly thicker than in the basin
	differenceFromExpectation <- sum(column$thickness) - sum(sedimentation)
	tolerance <- 0.01
	if (differenceFromExpectation>tolerance & differenceFromExpectation<=channelDepth) {
		cat("Column thickness is slightly thicker than what is shown in basin, likely owing to a channel at the base. The vertical coordinates of the column match what is shown in the basin, with 0.0 m being the base of the stratigraphy in the basin.\n")
	}

	column <- data.frame(modelTime=modelTime, thickness=column$thickness, stratPosition=column$stratPosition, facies=column$facies, elevation)
	
	class(column) <- "column"
	return(column)
}

#' @rdname column
#' @export

plot.column <- function(x, stratRange=c(floor(min(x$stratPosition)), ceiling(max(x$stratPosition))), axes=TRUE, ...) {
	# x is the output of the column() function	
	
	# Facies widths
	marineWidth <- 1.5
	carbonateShallowWidth <- 2.0
	carbonateDeepWidth <- 0.3
	floodplainWidth <- 1.5
	channelWidth <- 2
	sbChannelWidth <- 2.1
	paleosolWidth <- floodplainWidth
	unconformityWidth <- floodplainWidth
	plotWidth <- channelWidth
	
	# Facies colors
	marineColor <- 'tan'
	marineBorder <- 'tan3'
	carbonateColor <- 'cadetblue2'
	carbonateBorder <- 'cadetblue3'
	floodplainColor <- 'olivedrab3'
	floodplainBorder <- 'olivedrab4'
	channelColor <- 'yellow'
	sbChannelColor <- 'orangered'
	lagColor <- 'black'
	paleosolColor <- 'firebrick4'
	unconformityColor <- 'firebrick1'
	boxColor <- 'gray60'
	
	# Extract depositional units
	facies <- x$facies[-1]
	top <- x$stratPosition[-1]
	base <- x$stratPosition[-length(x$stratPosition)]
		
	# dev.new(height=7, width=3)
	plot(1, 1, type='n', xlim=c(0, plotWidth), ylim=stratRange, xlab='', ylab='Stratigraphic Position (m)', axes=FALSE, ...)
	if (axes == TRUE) {
		graphics::axis(2, las=1)
	}
	
	# Draw marine
	marine <- which(facies == 'marine')
	if (length(marine) > 0) {
		graphics::rect(rep(0, length(marine)), base[marine], rep(marineWidth, length(marine)), top[marine], col=marineColor, border=marineBorder)
	}
	
	# Draw carbonate
	carbonate <- which(facies == 'carbonate')
	if (length(carbonate) > 0) {
		cat("carbonates found\n")
		depth <- -x$elevation[carbonate]
		maxDepth <- max(depth)
		minDepth <- 0
		carbonateWidth <- rep(minDepth, length(carbonate))
		submerged <- depth > 0
		carbonateWidth[submerged] <- depth[submerged] * (carbonateShallowWidth - carbonateDeepWidth) / (maxDepth - minDepth) + carbonateDeepWidth	
		graphics::rect(rep(0, length(carbonate)), base[carbonate], carbonateWidth, top[carbonate], col=carbonateColor, border=carbonateBorder)
		graphics::rect(0, 0, carbonateDeepWidth, max(x$stratPosition), col="white", border="white", lwd=1.5)
		hiatus <- depth < 0
		graphics::segments(rep(0, length(carbonate[hiatus])), base[carbonate[hiatus]], carbonateDeepWidth, base[carbonate[hiatus]])
	}
	
	# Draw floodplain
	floodplain <- which(facies == 'floodplain')
	if (length(floodplain) > 0) {
		graphics::rect(rep(0, length(floodplain)), base[floodplain], rep(floodplainWidth, length(floodplain)), top[floodplain], col=floodplainColor, border=floodplainBorder)
	}
	
	# Draw channels
	channel <- which(facies == 'channel')
	if (length(channel) > 0) {
		# Channel bodies
		graphics::rect(rep(0, length(channel)), base[channel], rep(channelWidth, length(channel)), top[channel], col=channelColor, border=boxColor)

		# Channel lags
		graphics::segments(0, base[channel], channelWidth, base[channel], col=lagColor)
	}
	
	# Draw channels
	sbChannel <- which(facies == 'sequenceBoundingChannel')
	if (length(sbChannel) > 0) {
		# Channel bodies
		graphics::rect(rep(0, length(sbChannel)), base[sbChannel], rep(sbChannelWidth, length(sbChannel)), top[sbChannel], col=sbChannelColor, border=boxColor)

		# Channel lags
		graphics::segments(0, base[sbChannel], sbChannelWidth, base[sbChannel], col=lagColor)
	}
	
	# Draw unconformities 
	unconformity <- which(facies == 'unconformity')
	if (length(unconformity) > 0 ) {
		graphics::segments(0, base[unconformity], unconformityWidth, base[unconformity], col=unconformityColor, lwd=3.0)
	}
	
    # Draw paleosols
	paleosol <- which(facies == 'paleosol')
	if (length(paleosol) > 0) {
		graphics::segments(0, base[paleosol], paleosolWidth, base[paleosol], col=paleosolColor, lwd=3.0)
	}
}

#' 
#' @rdname column
#' @export

summary.column <- function(object, ...) {
	cat("duration:  ", max(object$modelTime), "m.y.\n")
	cat("thickness: ", sum(object$thickness), "m.y.\n")
}
