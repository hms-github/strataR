# Marine Module

# This offers similar functionality to the code in speciesModule.R and occurrencesModule.R, but entirely for marine species. Most of the functions have parallel names, but beginning with marine (for example marineAppendEcologicalParameters() vs. appendEcologicalParameters).

# Where no parallel function is created, the equivalent code in occurrencesModule.R should work.

# ---------------------------------------------------------------------------------------
## Marine species functions (paralleling those in speciesModule.R)

# Creating species is a two-step process. First, let species originate and go extinct through a random branching model (randomBranchTimeVarying). Second append the ecological  characteristics of those species.

# Use randomBranchTimeVarying() from speciesModule.R to simulate the marine species, as there is nothing environmentally specific in that code.

# appendEcologicalParameters(): Add preservation parameters to a species file

# id, ancestor: id numbers of a species and its ancestor; supplied in the data file
# origination, extinction: times of origination and speciation; supplied in the data file generated in the evolution module (e.g., by the randomBranchTimeVarying() function)

# minPD, maxPD: minimum and maximum preferred depth; values of minimum preferred depth less than 0 (above sea level) will cause distributions that are truncated at the shore, and this is generally necessary to keep diversity from dropping towards the coast
# meanDT, sdDT: mean and standard deviation of depth tolerance
# meanLogPA, sdLogPA: mean and standard deviation of the natural log of peak abundance, with peak abundance expressed as a percentage; Natural logs are used because PA is strongly right-skewed
# maxPA: maximum value of peak abundance, expressed as a percentage, to prevent excessively abundant species

# if meanLogPA and sdLogPA are set non-optimally, generating values frequently larger than maxPA, the while loop can run for a very long time, as it continually tries to reject PA values that are too large

marineAppendEcologicalParameters <- function(data, minPD, maxPD, meanDT, sdDT, meanLogPA, sdLogPA, maxPA) {
	
	# extract columns
	id          <- data$id
	ancestor    <- data$ancestor
	origination <- data$origination
	extinction  <- data$extinction
	numSpecies = length(ancestor)
	
	PD  <- runif(numSpecies, min=minPD, max=maxPD)
	DT  <- rnorm(numSpecies, mean=meanDT, sd=sdDT)
	PA <- exp(rnorm(numSpecies, mean=meanLogPA, sd=sdLogPA))
	
	# prevent excessively large values of PA
	while(max(PA) > maxPA) {
		numTooLarge <- sum(PA > maxPA)
		PA[PA > maxPA] <- exp(rnorm(numTooLarge, mean=meanLogPA, sd=sdLogPA))
	}
	
	# round all values
	PD  <- round(PD, 1)
	DT  <- round(DT, 1)
	PA  <- round(PA, 1)
	
	# create data frame
	species <- as.data.frame(cbind(id, ancestor, origination, extinction, PD, DT, PA))
	species
}



marineSpeciesCharacteristicsPlot <- function(species) {
	oldpar <- par(mfrow=c(3, 2), mar=c(5, 4, 0, 2) + 0.1)
	hist(species$origination, xlab='Time of origination (m.y.)', breaks=50, las=1, col='gray', main='')
	hist(species$extinction, xlab='Time of extinction (m.y.)', breaks=50, las=1, col='gray', main='')
	hist(species$PD, xlab='Preferred environment (m)', breaks=50, las=1, col='gray', main='')
	hist(species$DT, xlab='Depth tolerance (m)', breaks=50, las=1, col='gray', main='')
	hist(species$PA, xlab='Peak Abundance (as %)', breaks=50, las=1, col='gray', main='')
	par(oldpar)
}



# ---------------------------------------------------------------------------------------
## Marine occurrence functions (paralleling those in occurrencesModule.R)

# Probability of collection for a species at a given depth
# Used by marineGenerateOccurrences() and the diversity vs. water depth plots
marineProbCollection <- function(PD, DT, PA, waterDepth) {
	# PA is on a percent scale: 0-100
	probCollection <- (PA * exp(-((waterDepth - PD)^2)/(2 * DT^2))) / 100
	if (waterDepth <= 0) {           # marine (nonmarine would have a negative waterDepth)
		probCollection <- 0
	}
	probCollection
}



# Generate occurrences from stratigraphic history and species data
marineGenerateOccurrences <- function(stratColumn, species, sampleSpacing=0.5) {
	# NOTE: because columns are in terms of waterDepth, anything that is below sea level will have a negative value, so must be converted to a positive value to convert to water depth
	attach(stratColumn)  # modelTime, elevation, stratPosition, facies
	                                  # unused: thickness
	attach(species)  # id, ancestor, origination, extinction, PD, DT, PA
	
	numSpecies <- length(id)
	
	# vectors to hold results; oversize, then cull at end if necessary
	# age and position are equivalent to modelTime and waterDepth, will be renamed when
	#   data frame is created, but not here to avoid collision with attached strat data
	maxLength <- 100000
	storage <- vector(mode='numeric', length=maxLength)
	age       <- storage
	position  <- storage
	speciesId <- storage
	
	collectionHorizons <- seq(0, max(stratPosition), sampleSpacing)
	occurrenceNum <- 1
	
	for (horizon in collectionHorizons) {
		horizonTime <- modelTimeForHorizon(stratColumn, horizon, oldestAge=FALSE)
		horizonWaterDepth <- -elevation[horizonTime == modelTime]  # change sign to convert elevation to water depth
		for (i in 1:numSpecies) {
			if (horizonTime >= origination[i] && horizonTime <= extinction[i]) { # extant
				probCollection <- marineProbCollection(PD=PD[i], DT=DT[i], PA=PA[i], waterDepth=horizonWaterDepth)						
				random <- runif(1)
				if (random <= probCollection) {	                         # found
					age[occurrenceNum] <- horizonTime
					position[occurrenceNum] <- horizon
					speciesId[occurrenceNum] <- id[i]
					occurrenceNum <- occurrenceNum + 1
					if (occurrenceNum%%maxLength == 0) {
						age <- c(age, storage)
						position <- c(position, storage)
						speciesId <- c(speciesId, storage)
						message(paste("Resizing storage to hold up to", length(speciesId), "occurrences. If this happens multiple times, this may be a sign that the run needs to be terminated."))
					}
				}
			}
		}
	}	
	
	detach(stratColumn)
	detach(species)
	occurrences <- as.data.frame(cbind(modelTime=age, stratPosition=position, speciesId))
	# trim any rows with no data
	occurrences <- occurrences[-which(occurrences$speciesId == 0), ]
	message(paste("Generated", nrow(occurrences), "occurrences"))
	occurrences
}



# Plot all occurrences, along with a stratigraphic column and the history of waterDepth
marineOccurrenceColumnWaterDepthPlot <- function(occurrences, stratColumn, species, peakExtTimeMy=NA, extDurationMy=NA, orderedByLad=TRUE) {
	opar <- par(no.readonly = TRUE)
	leftDivider <- 0.20   # separates ranges from waterDepth plot
	rightDivider <- 0.80  # separates waterDepth plot from strat column
	
	# Left panel: stratigraphic column
	par(fig=c(0, leftDivider, 0, 1), mar=c(5, 4, 4, 0) + 0.1)
	stratColumnPlot(stratColumn)
	
	# Middle panel: range chart
	par(fig=c(leftDivider, rightDivider, 0, 1), mar=c(5, 0, 4, 0) + 0.1, new=TRUE)
	occurrencePlot(occurrences=occurrences, stratColumn=stratColumn, species=species, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy, orderedByLad=orderedByLad, yAxisLabels=FALSE)
	
	# Right panel: water depth plot
	par(fig=c(rightDivider, 1, 0, 1), mar=c(5, 0, 4, 2) + 0.1, new=TRUE)
	waterDepthHistoryPlot(stratColumn, yAxisLabels=FALSE)
	
	par(opar)
}



# Plot diversity vs. water depth, used as a diagnostic plot to see if diversity covaries with water depth (it generally shouldn't)
# In this function, a species is considered to occur if its probability of collection exceeds a specified threshold (pCrit)
# waterDepthDiversityResamplePlot() is the preferred way to do this, as it adds a confidence interval, but it is also slower
marineWaterDepthDiversityPlot <- function(species, minWaterDepth=0, maxWaterDepth=150, pCrit=0.1) {
	attach(species)
	numSpecies <- nrow(species)
	
	waterDepth <- seq(minWaterDepth, maxWaterDepth, 1)
	diversity <- rep(0, length(waterDepth))
	
	for (depth in 1:length(waterDepth)) {
		waterDepthDiversity <- 0
		for (i in 1:numSpecies) {
			marineProbCollection <- marineProbCollection(PD=PD[i], DT=DT[i], PA=PA[i], waterDepth=waterDepth[depth])
			if (marineProbCollection >= pCrit) {	                         # found
				waterDepthDiversity <- waterDepthDiversity + 1
			}
		}
		diversity[depth] <- waterDepthDiversity
	}
	
	detach(species)
	plot(waterDepth, diversity, ylim=c(0, max(diversity)), xlab='Water Depth (m)', ylab='Richness (number of species)', type='l', lwd=2, las=1)
}



# Diagnostic plot to evaluate whether diversity changes with water depth (it generally shouldn't)
# In this function, a species is considered to occur if a random number is less than its probability of collection.
# This is done over a number of trials to place a confidence interval on the observed diversity at any water depth.
# Function is slow, owing to the three nested loops.
marineWaterDepthDiversityResamplePlot <- function(species, minWaterDepth=0, maxWaterDepth=150, time=0, numTrials=100) {
	selectedSpecies <- species[which(species$origination <= time & species$extinction >= time), ]
	numSpecies <- nrow(selectedSpecies)
	attach(selectedSpecies)
	
	waterDepth <- seq(minWaterDepth, maxWaterDepth, 1)
	diversity50 <- rep(0, length(waterDepth)) # mean
	diversity25 <- rep(0, length(waterDepth)) # 25% quantile
	diversity75 <- rep(0, length(waterDepth)) # 75% quantile
	
	for (depth in 1:length(waterDepth)) {
		trialDiversity <- rep(0, numTrials)
		for (trial in 1:numTrials) {
			iterationDiversity <- 0
			for (i in 1:numSpecies) {
				random <- runif(1)
				marineProbCollection <- marineProbCollection(PD=PD[i], DT=DT[i], PA=PA[i], waterDepth=waterDepth[depth])
				if (random <= marineProbCollection) {	                         # found
					iterationDiversity <- iterationDiversity + 1
				}
			}
			trialDiversity[trial] <- iterationDiversity
		}
		diversity50[depth] <- quantile(trialDiversity, 0.50)
		diversity25[depth] <- quantile(trialDiversity, 0.25)
		diversity75[depth] <- quantile(trialDiversity, 0.75)
	}

	detach(selectedSpecies)
	yLimits <- c(0, max(diversity75))
	plot(waterDepth, diversity50, ylim=yLimits, xlab='Water Depth (m)', ylab='Richness (number of species)', type='n', lwd=2, las=1)

	polyX <- c(waterDepth, rev(waterDepth))
	polyY <- c(diversity25, rev(diversity75))
	polygon(polyX, polyY, col='lightgray', border='lightgray')
	
	points(waterDepth, diversity50, type='l', lwd=2, col='black')
}

# ---------------------------------------------------------------------------------------
## SPECIES PLOTS

marineSpeciesCharacteristicsPlot <- function(species) {
	oldpar <- par(mfrow=c(3, 2), mar=c(5, 4, 0, 2) + 0.1)
	hist(species$origination, xlab='Time of origination (m.y.)', breaks=50, las=1, col='gray', main='')
	hist(species$extinction, xlab='Time of extinction (m.y.)', breaks=50, las=1, col='gray', main='')
	hist(species$PD, xlab='Preferred depth (m)', breaks=50, las=1, col='gray', main='')
	hist(species$DT, xlab='Depth tolerance (m)', breaks=50, las=1, col='gray', main='')
	hist(species$PA, xlab='Peak Abundance (as %)', breaks=50, las=1, col='gray', main='')
	par(oldpar)
}



# ---------------------------------------------------------------------------------------
# Sedflux utilities

# Converts a sedflux column to one that is compatible with the nonmarine model
convertSedfluxColumn <- function(x, timeMultiplier) {
	colnames(x) <- c('modelTime', 'elevation', 'stratPosition')
	x$elevation <- - x$elevation # flip sign so negative elevation is now positive water depth
	x$modelTime <- x$modelTime / 1000000 # convert time from years to millions of years
	x$modelTime <- x$modelTime * timeMultiplier
	return(x)
}

# Function to remove eroded strata from a sedflux history
# Necessary to prevent loops in a stratigraphic position vs. water depth plot, and to
#  avoid sampling horizons present in the history file that are not actually preserved.

sedfluxCullErodedLayers <- function(fileName) {
	data <- read.table(fileName, header=TRUE, sep=',')
	message(paste('Read file', fileName))
	message('Now culling eroded layers; this may take a couple of minutes\n')
	i <- length(data$stratPosition)
	while(i>1) {
		if (data$stratPosition[i] > data$stratPosition[i-1]) {			# deposition
			i <- i-1  # go to next layer
		} else if (data$stratPosition[i] < data$stratPosition[i-1]) {	# erosion
			# need to check if the entire layer is eroded or if it is partly preserved
			if (data$stratPosition[i] < data$stratPosition[i-2]) {	# entire i-1 layer eroded
				erodedLayer <- i-1
				data <- data[-erodedLayer, ]
				i <- i-1  # update my layer index
			} else {					# partial erosion of i-1 layer
				amountEroded <- abs(data$stratPosition[i] - data$stratPosition[i-1])
				data$stratPosition[i-1] <- data$stratPosition[i-1] - amountEroded
			}
		} else {													# no deposition or erosion
			# remove the current horizon
			data <- data[-i, ]
			i <- i-1  # to to next layer
		}
	}
	filePrefix <- paste(head(unlist(strsplit(fileName, split='.', fixed=TRUE)), -1), collapse='')
	writtenFileName <- paste(filePrefix, '-culled.csv', sep='')
	write.table(data, file=writtenFileName, sep=',', row.names=FALSE)
	message(paste('Culled history saved as', writtenFileName))
}

