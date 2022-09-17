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
	
	PD  <- stats::runif(numSpecies, min=minPD, max=maxPD)
	DT  <- stats::rnorm(numSpecies, mean=meanDT, sd=sdDT)
	PA <- exp(stats::rnorm(numSpecies, mean=meanLogPA, sd=sdLogPA))
	
	# prevent excessively large values of PA
	while(max(PA) > maxPA) {
		numTooLarge <- sum(PA > maxPA)
		PA[PA > maxPA] <- exp(stats::rnorm(numTooLarge, mean=meanLogPA, sd=sdLogPA))
	}
	
	# round all values
	PD  <- round(PD, 1)
	DT  <- round(DT, 1)
	PA  <- round(PA, 1)
	
	# prevent negative DT (assigned an arbitrarily small positive value)
	if (min(DT) <= 0) {
		DT[DT<=0] <- 0.1
		message("WARNING: Negative values of DT generated; replaced with 0.1")
	}
	
	# create data frame
	species <- as.data.frame(cbind(id, ancestor, origination, extinction, PD, DT, PA))
	species
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
	data <- utils::read.table(fileName, header=TRUE, sep=',')
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
	filePrefix <- paste(utils::head(unlist(strsplit(fileName, split='.', fixed=TRUE)), -1), collapse='')
	writtenFileName <- paste(filePrefix, '-culled.csv', sep='')
	utils::write.table(data, file=writtenFileName, sep=',', row.names=FALSE)
	message(paste('Culled history saved as', writtenFileName))
}

