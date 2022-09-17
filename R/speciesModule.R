# Species Module

# ---------------------------------------------------------------------------------------
## SIMULATE SPECIES EVOLUTION AND ECOLOGY

# Creating species is a two-step process. First, let species originate and go extinct through a random branching model (randomBranchTimeVarying). Second append the ecological  characteristics of those species.

# randomBranchTimeVarying(): Generate species using a random branching model that allows for a mass extinction

# all time units are in millions of years

# timeStep should be set to a fairly short value (e.g., 0.0001 m.y.)
# extRate is the extinction rate per million years, set by default to the Phanerozoic average of 0.25 per million years, from Raup
# startingSpecies is the number of species at the beginning of the random branching model
# durationMy is the total time span of the simulation

# By default, the random-branching model is time-homogeneous. By changing peakExtFactor, peakExtDurationMy, and peakExtTimeMy, one can add a period of elevated extinction rates. 
# peakExtFactor is a multiplier, describing how many much greater the peak extinction rates during the mass extinction are above background (e.g., 25x). The default value of 1 stipulates that there is no mass extinction.
# peakExtTimeMy is the moment in time where extinction rates achieve the peakExtFactor.
# peakExtDurationMy is the total time span of the extinction. During this time, extinction rate will be rise (or fall) linearly between the background and the peak extinction rate. 	
	
randomBranchTimeVarying <- function(timeStep, extRate=0.25, startingSpecies=200, durationMy=10, peakExtFactor=1, peakExtDurationMy=1, peakExtTimeMy=2) {
	
	# all calculations are in millions of years
	distantFuture <- durationMy + 1    # run 1 million years longer

	# per time step probability of extinction and origination
	extinctionProb <- extRate * timeStep
	speciationProb <- extinctionProb
	peakExtRate <- extRate * peakExtFactor
	
	modelTime <- seq(0, durationMy, timeStep)
	currentExtRate <- dnorm(modelTime, mean=peakExtTimeMy, sd=peakExtDurationMy/6)   # duration = 6 s.d.
	currentExtRate <- currentExtRate/max(currentExtRate)*(peakExtRate-extRate) + extRate
	
	# vectors to hold times of origination and extinction, ancestors
	origination <- rep(0, startingSpecies)
	extinction <- rep(distantFuture, startingSpecies)
	ancestor <- rep(0, startingSpecies)
	id <- 1:startingSpecies
	
	nextSpecies <- startingSpecies + 1
	numTimeSteps <- length(modelTime)
	for (tStep in 1:numTimeSteps) {
		timeNow <- modelTime[tStep]
		extRateNow <- currentExtRate[tStep] * timeStep
		currentSpecies <- nextSpecies - 1
		for (i in 1:currentSpecies) {
			if (extinction[i] > timeNow) {					# still extant
				if (runif(1) < extRateNow) {				# went extinct
					extinction[i] <- timeNow
				}
				if (runif(1) < speciationProb) {			# speciated
					origination[nextSpecies] <- timeNow
					extinction[nextSpecies]  <- distantFuture
					ancestor[nextSpecies]    <- i
					id[nextSpecies] <- nextSpecies
					nextSpecies <- nextSpecies + 1
				}
			}
		}
	}
	
	species <- as.data.frame(cbind(id, ancestor, origination, extinction))
	message(paste("Generated", nrow(species), "species"))
	species
}

# appendEcologicalParameters(): Add preservation parameters to a species file

# id, ancestor: id numbers of a species and its ancestor; supplied in the data file
# origination, extinction: times of origination and speciation; supplied in the data file generated in the evolution module (e.g., by the randomBranchTimeVarying() function)

# minPE, maxPE: minimum and maximum preferred elevation; values of minPE less than 0 (sea level) will cause distributions that are truncated at the shore, and this is generally necessary to keep diversity from dropping towards the coast
# meanET, sdET: mean and standard deviation of elevation tolerance
# meanLogPA, sdLogPA: mean and standard deviation of the natural log of peak abundance, with peak abundance expressed as a percentage; Natural logs are used because PA is strongly right-skewed
# maxPA: maximum value of peak abundance, expressed as a percentage, to prevent excessively abundant species
# Affinity is the preference of species for being preserved in channels in floodplain, and it ranges from -1 to 1. -1 is a species preserved only in channels, +1 is a species preserved only on floodplains, and 0 is a species is preserved equally in both
# The defaults for affinity, minAff and maxAff, span this entire range

# if meanLogPA and sdLogPA are set non-optimally, generating values frequently larger than maxPA, the while loop can run for a very long time, as it continually tries to reject PA values that are too large

appendEcologicalParameters <- function(data, minPE, maxPE, meanET, sdET, meanLogPA, sdLogPA, maxPA, minAff=-1, maxAff=1) {
	
	# extract columns
	id          <- data$id
	ancestor    <- data$ancestor
	origination <- data$origination
	extinction  <- data$extinction
	numSpecies = length(ancestor)
	
	PE  <- runif(numSpecies, min=minPE, max=maxPE)
	ET  <- rnorm(numSpecies, mean=meanET, sd=sdET)
	PA <- exp(rnorm(numSpecies, mean=meanLogPA, sd=sdLogPA))
	Aff <- runif(numSpecies, min=minAff, max=maxAff)
	
	# prevent excessively large values of PA
	while(max(PA) > maxPA) {
		numTooLarge <- sum(PA > maxPA)
		PA[PA > maxPA] <- exp(rnorm(numTooLarge, mean=meanLogPA, sd=sdLogPA))
	}
	
	# round all values
	PE  <- round(PE, 1)
	ET  <- round(ET, 1)
	PA  <- round(PA, 1)
	Aff <- round(Aff, 2)
	
	# create data frame
	species <- as.data.frame(cbind(id, ancestor, origination, extinction, PE, ET, PA, Aff))
	species
}




# ---------------------------------------------------------------------------------------
## SPECIES PLOTS

speciesCharacteristicsPlot <- function(species) {
	oldpar <- par(mfrow=c(3, 2), mar=c(5, 4, 0, 2) + 0.1)
	hist(species$origination, xlab='Time of origination (m.y.)', breaks=50, las=1, col='gray', main='')
	hist(species$extinction, xlab='Time of extinction (m.y.)', breaks=50, las=1, col='gray', main='')
	hist(species$PE, xlab='Preferred environment (m)', breaks=50, las=1, col='gray', main='')
	hist(species$ET, xlab='Environmental tolerance (m)', breaks=50, las=1, col='gray', main='')
	hist(species$PA, xlab='Peak Abundance (as %)', breaks=50, las=1, col='gray', main='')
	hist(species$Aff, xlab='Affinity', breaks=50, las=1, col='gray', main='')
	par(oldpar)
}



