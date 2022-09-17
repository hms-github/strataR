# Sedimentation Module

# ---------------------------------------------------------------------------------------
## SIMULATE STRATIGRAPHIC COLUMNS

# Generate a stratigraphic column from a basin produced by fillBasin() in balancedFillModule
generateColumnFromBasin <- function(basin, locationKm, setting=c('valley', 'interfluve'), pChannel, channelDepth) {
	# basin: output of fillBasin()
	# locationKm: position in km from left edge of basin; nearest location will be used
	# pChannel: probability of a channel being present (ranges from 0.0 to 1.0)
	# channelDepth: in m
	
	setting <- match.arg(setting)
	if (setting != 'valley' && setting != 'interfluve') {
	     warning("setting must be 'valley' or 'interfluve'")
	}
	
	elevation <- elevationAtLocation(basin, locationKm, setting)             # m
	sedimentation <- sedimentationAtLocation(basin, locationKm, setting)     # m
	timeStep <- basin$parameters$timeStep                                    # m.y.
	
	# Create vectors to hold modelTime, thickness, stratPosition, and facies:
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
			} else if (runif(1) > pChannel) {				# floodplain
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
		message(paste("Column thickness is slightly thicker than what is shown in basin, likely owing to a channel at the base. The vertical coordinates of the column match what is shown in the basin, with 0.0 m being the base of the stratigraphy in the basin."))
	}

	column <- data.frame(modelTime=modelTime, thickness=column$thickness, stratPosition=column$stratPosition, facies=column$facies, elevation)
	return(column)
}

# ---------------------------------------------------------------------------------------

# Generate a stratigraphic column not tied to fillBasin()
generateColumn <- function(totalTime, timeStep, elevation, pChannel, channelDepth, subsidenceRate) {
	# totalTime: duration of simulation, in years
	# timeStep: length of a single time step, in years
	# elevation: vector of elevations at location, in m
	# pChannel: probability of a channel as a single value, in 0.0 – 1.0 range
	# channelDepth: single value, in m
	# subsidenceRate: vector with rates in m/yr
	
	# Create vectors to hold modelTime, thickness, stratPosition, and facies:
	modelTime     <- seq(0, totalTime, timeStep)
	thickness     <- rep(0,  length(modelTime))
	stratPosition <- rep(0,  length(modelTime))
	facies        <- rep('', length(modelTime))
	
	# Check that all vectors have lengths equal to the number of time steps
	numTimeSteps <- totalTime/timeStep + 1
	error <- FALSE
	if (length(elevation) != numTimeSteps) {
	  message('length of the elevation vector does not match the number of time steps')
	  error <- TRUE
	}
	if (length(subsidenceRate) != numTimeSteps) {
	  message('length of the subsidenceRate vector does not match the number of time steps')
	  error <- TRUE
	}
	if (length(modelTime) != numTimeSteps) {
	  message('length of the modelTime vector does not match the number of time steps')
	  error <- TRUE
	}
	if (length(thickness) != numTimeSteps) {
	  message('length of the thickness vector does not match the number of time steps')
	  error <- TRUE
	}
	if (length(stratPosition) != numTimeSteps) {
	  message('length of the stratPosition vector does not match the number of time steps')
	  error <- TRUE
	}
	if (length(facies) != numTimeSteps) {
	  message('length of the facies vector does not match the number of time steps')
	  error <- TRUE
	}
	
	if (error) {
		stop('generateColumn() aborted', call.=FALSE)
	}

	## Build column
	netIncision <- 0
	for (i in 2:length(modelTime)) { # first time step is base of column
		#calculate accommodation space
		elevationIncrease <- elevation[i] - elevation[i-1]
		subsidence <- subsidenceRate[i] * timeStep
		accommodation <- subsidence + elevationIncrease
		
		if (runif(1) > pChannel) {    # No channel
			if (accommodation > 0) {  #    accommodation, deposit floodplain
				thickness[i] <- accommodation
				stratPosition[i] <- thickness[i] + stratPosition[i-1]
				facies[i] <- 'floodplain'
				netIncision <- 0
			} else {                  #    no accommodation, subaerial exposure
				thickness[i] <- 0
				stratPosition[i] <- stratPosition[i-1]
				facies[i] <- 'paleosol'
			}  
		} else {                      # Channel is present
			amountToBeEroded <- channelDepth - accommodation
			if (amountToBeEroded <= 0) {  # accommodation, fill space with channel   
				thickness[i] <- accommodation
				stratPosition[i] <- thickness[i] + stratPosition[i-1]
				facies[i] <- 'channel'
			} else {                      # some erosion needed
				# erode to base of channel
				columnTop <- stratPosition[i-1]
				depthOfErosion <- columnTop - amountToBeEroded
				horizonsToBeRemoved <- which(stratPosition > depthOfErosion)
				horizonsToCompletelyRemove <- horizonsToBeRemoved[-1]

				# completely erode all but lowest horizon      
				thickness[horizonsToCompletelyRemove] <- 0 
				stratPosition[horizonsToCompletelyRemove] <- depthOfErosion
				facies[horizonsToCompletelyRemove] <- 'channel'

				# partly erode the lowest horizon
				remainingErosion <- stratPosition[horizonsToBeRemoved[1]] - depthOfErosion    
				thickness[horizonsToBeRemoved[1]] <- thickness[horizonsToBeRemoved[1]] - remainingErosion
				stratPosition[horizonsToBeRemoved[1]] <- depthOfErosion
				netIncision <- netIncision + amountToBeEroded

				# deposit the channel
				thickness[i] <- channelDepth                                         
				stratPosition[i] <- thickness[i] + stratPosition[i-1]
				if (netIncision <= channelDepth) {
					facies[i] <- 'channel'
				} else {
					facies[i] <- 'sequenceBoundingChannel'
				}
				
			}
		}
	}
	
	# Corrected thicknesses for when channel occurs near base of column (therefore 
	#   creating negative thicknesses near base of column)
	minThickness <- min(thickness)
	stratPosition <- stratPosition + abs(minThickness)
	thickness <- c(0, diff(stratPosition))
	
	column <- data.frame(modelTime, thickness, stratPosition, facies, elevation)
	return(column)
}




# ---------------------------------------------------------------------------------------
## STRATIGRAPHIC COLUMN PLOTS

# Plot a stratigraphic column
stratColumnPlot <- function(column, stratRange=c(floor(min(column$stratPosition)), ceiling(max(column$stratPosition))), axes=TRUE, ...) {
	# column is the output of the generateColumn() function	
	
	# Facies widths
	marineWidth <- 1.5
	floodplainWidth <- 1.5
	channelWidth <- 2
	sbChannelWidth <- 2.1
	paleosolWidth <- floodplainWidth
	unconformityWidth <- floodplainWidth
	plotWidth <- channelWidth
	
	# Facies colors
	marineColor <- 'tan'
	marineBorder <- 'tan3'
	floodplainColor <- 'olivedrab3'
	floodplainBorder <- 'olivedrab4'
	channelColor <- 'yellow'
	sbChannelColor <- 'orangered'
	lagColor <- 'black'
	paleosolColor <- 'firebrick4'
	unconformityColor <- 'firebrick1'
	boxColor <- 'gray60'
	
	# Extract depositional units
	facies <- column$facies[-1]
	top <- column$stratPosition[-1]
	base <- column$stratPosition[-length(column$stratPosition)]
		
	# dev.new(height=7, width=3)
	plot(1, 1, type='n', xlim=c(0, plotWidth), ylim=stratRange, xlab='', ylab='Stratigraphic Position (m)', axes=FALSE, ...)
	if (axes == TRUE) {
		axis(2, las=1)
	}
	
	# Draw marine
	marine <- which(facies == 'marine')
	if (length(marine) > 0) {
		rect(rep(0, length(marine)), base[marine], rep(marineWidth, length(marine)), top[marine], col=marineColor, border=marineBorder)
	}
	
	# Draw floodplain
	floodplain <- which(facies == 'floodplain')
	if (length(floodplain) > 0) {
		rect(rep(0, length(floodplain)), base[floodplain], rep(floodplainWidth, length(floodplain)), top[floodplain], col=floodplainColor, border=floodplainBorder)
	}
	
	# Draw channels
	channel <- which(facies == 'channel')
	if (length(channel) > 0) {
		# Channel bodies
		rect(rep(0, length(channel)), base[channel], rep(channelWidth, length(channel)), top[channel], col=channelColor, border=boxColor)

		# Channel lags
		segments(0, base[channel], channelWidth, base[channel], col=lagColor)
	}
	
	# Draw channels
	sbChannel <- which(facies == 'sequenceBoundingChannel')
	if (length(sbChannel) > 0) {
		# Channel bodies
		rect(rep(0, length(sbChannel)), base[sbChannel], rep(sbChannelWidth, length(sbChannel)), top[sbChannel], col=sbChannelColor, border=boxColor)

		# Channel lags
		segments(0, base[sbChannel], sbChannelWidth, base[sbChannel], col=lagColor)
	}
	
	# Draw unconformities 
	unconformity <- which(facies == 'unconformity')
	if (length(unconformity) > 0 ) {
		segments(0, base[unconformity], unconformityWidth, base[unconformity], col=unconformityColor, lwd=3.0)
	}
	
    # Draw paleosols
	paleosol <- which(facies == 'paleosol')
	if (length(paleosol) > 0) {
		segments(0, base[paleosol], paleosolWidth, base[paleosol], col=paleosolColor, lwd=3.0)
	}
}

# Plot the sedimentation rates, binned over a specified interval (in meters) of the column
sedRateHistoryPlot <- function(stratColumn, stratBin=10, yAxisLabels=TRUE, ...) {
	rates <- sedRateAveraged(stratColumn, stratBin=stratBin)
	ylab <- c('Stratigraphic Position (m)')
	if (yAxisLabels == FALSE) {
		ylab = c('')
	}
	yLimits <- c(min(stratColumn$stratPosition), max(stratColumn$stratPosition))
	plot(rates$sedRate / 1000, rates$stratPosition, las=1, xlab='Sed. Rate (m/kyr)', ylab=ylab, ylim=yLimits, type='o', col='brown', pch=16, axes=yAxisLabels, frame.plot=yAxisLabels, ...)
	if (yAxisLabels == FALSE) {
		axis(1)
	}
}

# Plot the elevation recorded in a stratigraphic column
elevationHistoryPlot <- function(stratColumn, yAxisLabels=TRUE, ...) {
	ylab <- c('Stratigraphic Position (m)')
	if (yAxisLabels == FALSE) {
		ylab = c('')
	}
	yLimits <- c(min(stratColumn$stratPosition), max(stratColumn$stratPosition))
	plot(stratColumn$elevation, stratColumn$stratPosition, type='n', xlab='Elevation (m)', ylab=ylab, ylim=yLimits, las=1, axes=yAxisLabels, frame.plot=yAxisLabels, ...)
	if (yAxisLabels == FALSE) {
		axis(1)
	}
	abline(v=0, col='gray', lty='dotted')
	points(stratColumn$elevation, stratColumn$stratPosition, type='l', lwd=1.5, col='chartreuse4')
}




# ---------------------------------------------------------------------------------------
## UTILITY FUNCTIONS

# Find the modelTime of a horizon in a stratigraphic column plot by clicking on it
horizonAge <- function(stratColumn) {
	# stratColumn should have at least 3 columns: elevation, stratPosition, modelTime
	plot(stratColumn$elevation, stratColumn$stratPosition, type='l')
	point <- locator(1)
	selectedColumnPosition <- point$y[1]
	diff <- abs(stratColumn$stratPosition - selectedColumnPosition)
	row <- which(diff == min(diff))
	modelTime <- stratColumn$modelTime[row]
	modelTime
}

# Find the modelTime of a specified horizon in a stratigraphic column
# Used by sedRateAveraged(), and by generateOccurrences() in occurrencesModule
modelTimeForHorizon <- function(stratColumn, horizon, oldestAge=TRUE) {
	diff <- abs(stratColumn$stratPosition - horizon)
	row <- 0
	if (oldestAge == TRUE) {
		row <- min(which(diff == min(diff)))
	} else {
		row <- max(which(diff == min(diff)))
	}
	modelTime <- stratColumn$modelTime[row]
	return(modelTime)
}

# Find the stratigraphic position of a specified model time in a stratigraphic column 
# Used by occurrencePlot() in occurrencesModule
stratPositionForAge <- function(modelTime, stratColumn) {
	stratPosition <- -9999
	if (modelTime > max(stratColumn$modelTime)) {
		stratPosition <- 9999
	} else if (modelTime < min(stratColumn$modelTime)) {
		stratPosition <- -9999
	} else {
		row <- min(which(stratColumn$modelTime >= modelTime))
		stratPosition <- stratColumn$stratPosition[row]
	}
	stratPosition
}


# Find the sedimentation rate through a column for equal-sized bins (stratBin)
# Used by sedRateHistoryPlot()
sedRateAveraged <- function(stratColumn, stratBin=10) {
	# units of stratBin are meters
	binBoundary <- seq(0, max(stratColumn$stratPosition), stratBin)
	binBoundary <- c(binBoundary, max(binBoundary) + stratBin)
	oldestAge <- rep(0, length(binBoundary))
	for (i in 1:length(binBoundary)) {
		oldestAge[i] <- modelTimeForHorizon(stratColumn=stratColumn, horizon=binBoundary[i], oldestAge=TRUE)
	}
	duration <- diff(oldestAge)        # m.y.
	sedRate <- stratBin / duration     # m / m.y.
	binBoundary <- binBoundary[-1]
	
	# remove last point
	lastPoint <- length(binBoundary)
	binBoundary <- binBoundary[-lastPoint]
	sedRate <- sedRate[-lastPoint]
	
	rates <- data.frame(cbind(stratPosition=binBoundary, sedRate=sedRate))
	return(rates)
}


# ---------------------------------------------------------------------------------------
# HELPER FUNCTIONS USED BY generateColumnFromBasin(), GENERALLY NOT TO BE CALLED DIRECTLY BY USERS

depositMarine <- function(column, timeStep) {
	column$thickness[timeStep] <- column$sedimentation[timeStep]
	column$stratPosition[timeStep] <- column$stratPosition[timeStep-1] + column$sedimentation[timeStep]
	column$facies[timeStep] <- 'marine'
	return(column)
}

depositFloodplain <- function(column, timeStep) {
	column$thickness[timeStep] <- column$sedimentation[timeStep]
	column$stratPosition[timeStep] <- column$stratPosition[timeStep-1] + column$sedimentation[timeStep]
	column$facies[timeStep] <- 'floodplain'
	return(column)
}

depositNonincisingChannel <- function(column, timeStep) {
	column$thickness[timeStep] <- column$sedimentation[timeStep]
	column$stratPosition[timeStep] <- column$stratPosition[timeStep-1] + column$sedimentation[timeStep]
	column$facies[timeStep] <- 'channel'
	return(column)
}

depositIncisingChannel <- function(column, timeStep, channelDepth) {
	column$thickness[timeStep] <- channelDepth
	column$stratPosition[timeStep] <- column$stratPosition[timeStep-1] + column$sedimentation[timeStep]
	column$facies[timeStep] <- 'channel'

	# now correct older deposits, owing to erosion
	remainingIncision <- channelDepth - column$sedimentation[timeStep]
	depthOfErosion <- column$stratPosition[timeStep] - channelDepth
	underlyingHorizons <- 1:(timeStep-1)
	horizonsToBeRemoved <- which(column$stratPosition[underlyingHorizons] >= depthOfErosion)				
	unconformityPresent <- 'unconformity' %in% column$facies[horizonsToBeRemoved]
	horizonsToCompletelyRemove <- horizonsToBeRemoved[-1]
	
	lowestHorizon <- 1
	if (lowestHorizon %in% horizonsToBeRemoved) { # eroding through base of basin
		column$thickness[timeStep] <- channelDepth
		column$thickness[horizonsToBeRemoved] <- 0
		column$stratPosition[timeStep] <- column$sedimentation[timeStep]
		column$stratPosition[horizonsToBeRemoved] <- column$sedimentation[timeStep] - channelDepth
		column$facies[horizonsToBeRemoved] <- 'channel'
		# correct the data itself (sedimentation) so that basin and column agree
		column$sedimentation[horizonsToBeRemoved] <- 0
	} else {  # this is the normal case
		# completely erode all but lowest horizon      
		column$thickness[horizonsToCompletelyRemove] <- 0 
		column$stratPosition[horizonsToCompletelyRemove] <- depthOfErosion
		column$facies[horizonsToCompletelyRemove] <- 'channel'
	
		# partly erode the lowest horizon
		remainingErosion <- column$stratPosition[horizonsToBeRemoved[1]] - depthOfErosion    
		column$thickness[horizonsToBeRemoved[1]] <- column$thickness[horizonsToBeRemoved[1]] - remainingErosion
		column$stratPosition[horizonsToBeRemoved[1]] <- depthOfErosion
		if (unconformityPresent && length(horizonsToBeRemoved)>=3) {
			column$facies[horizonsToBeRemoved[2]] <- 'unconformity'
			column$stratPosition[horizonsToBeRemoved[2]] <- depthOfErosion
		}
	}
	
	return(column)
}

nondepositionDuringBypass <- function(column, timeStep) {
	column$thickness[timeStep] <- 0
	column$stratPosition[timeStep] <- column$stratPosition[timeStep-1]
	column$facies[timeStep] <- 'paleosol'
	return(column)
}

erodeAtUnconformity <- function(column, timeStep) {
	column$thickness[timeStep] <- 0
	column$stratPosition[timeStep] <- column$stratPosition[timeStep-1] + column$sedimentation[timeStep]  # erosion
	column$facies[timeStep] <- 'unconformity'

	# now correct older deposits, owing to erosion
	remainingIncision <- -column$sedimentation[timeStep]
	depthOfErosion <- column$stratPosition[timeStep]
	underlyingHorizons <- 1:(timeStep-1)
	horizonsToBeRemoved <- which(column$stratPosition[underlyingHorizons] >= depthOfErosion)
	horizonsToCompletelyRemove <- horizonsToBeRemoved[-1]
	
	lowestHorizon <- 1
	if (lowestHorizon %in% horizonsToBeRemoved) { # eroding through base of basin
		column$thickness[timeStep] <- 0
		column$thickness[horizonsToBeRemoved] <- 0
		column$stratPosition[timeStep] <- 0
		column$stratPosition[horizonsToBeRemoved] <- 0
		column$facies[horizonsToBeRemoved] <- 'unconformity'
		column$facies[1] <- ''
		# correct the data itself (sedimentation) so that basin and column agree
		column$sedimentation[timeStep] <- 0
		column$sedimentation[horizonsToBeRemoved] <- 0
	} else {      # this is the normal case
		# completely erode all but lowest horizon      
		column$thickness[horizonsToCompletelyRemove] <- 0 
		column$stratPosition[horizonsToCompletelyRemove] <- depthOfErosion
		column$facies[horizonsToCompletelyRemove] <- 'unconformity'

		# partly erode the lowest horizon
		remainingErosion <- column$stratPosition[horizonsToBeRemoved[1]] - depthOfErosion    
		column$thickness[horizonsToBeRemoved[1]] <- column$thickness[horizonsToBeRemoved[1]] - remainingErosion
		column$stratPosition[horizonsToBeRemoved[1]] <- depthOfErosion
	}
	
	return(column)
}

# ---------------------------------------------------------------------------------------
# DEAD CODE
# Code that is apparently not called by any other function or not called by user
# Consider archiving

# Find the instantaneous sedimentation rate for all horizons in a stratigraphic column
# [DEAD]
sedimentationRate <- function(stratColumn) {	
	rate <- c(diff(stratColumn$stratPosition)/diff(stratColumn$modelTime), 0)
	rate
}

# [DEAD]
testEquals <- function(x, y) {
	if (identical(x, y)) {
		print('Values are equal')
	} else {
		warning(paste('WARNING: Test failed. Values are not equal.', x, ' vs. ', y))
	}
}

# [DEAD]
testEqualWithinTolerance <- function(x, y, tol=0.00000001) {
	message <- 'Values are equal'
	issueWarning <- FALSE
	n <- length(x)
	for (i in 1:n) {
		if (x[i] > y[i] + tol) {
			issueWarning <- TRUE
			message <- 'WARNING: values are not all equal'
		} else if (x[i] < y[i] - tol) {
			issueWarning <- TRUE
			message <- 'WARNING: values are not all equal'
		}
	}
	if (issueWarning) {
		warning('WARNING: values are not all equal')
	} else {
		print('Values are equal')
	}
}



