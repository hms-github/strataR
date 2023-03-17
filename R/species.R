#' @title Species Simulation
#' 
#' @description Create an object of class "species"
#'
#' @details Creates a species object, which records the ancestor of a species, its times of origination and extinction, and its ecological characteristics: preferred elevation, elevation tolerance, peak abundance, and channel/floodplain affinity. Evolution of species is treated as a random-branching process, which can be optionally be decreased during a mass extinction. A species object is designed for nonmarine species; for marine species, see marineSpecies.
#' 
#' @param timeStep the time step for simulating evolution, in m.y. Typically set to a  
#'   small value, such as 0.001 (1000 years).
#' @param extRate extinction rate, as a probability per m.y. Speciation rate is set to be  
#'   equal to the extinction rate, such that there is no forced change in diversity other  
#'   than what occurs through drift.
#' @param startingSpecies number of species at start of simulation.
#' @param durationMy duration of species simulation, in m.y.
#' @param peakExtFactor factor (multiplier) describing how much extinction rate will  
#'   increase at the peak of a mass extinction. Set to 1 to have no mass extinction.
#' @param peakExtDurationMy duration of mass extinction, in m.y. 
#' @param peakExtTimeMy time at which mass extinction hits its peak rate, in m.y.  
#'   Extinction rate will climb linearly from background rates to peak rates through the  
#'   first half of the extinction window, then decrease linearly back to background rates  
#'   in the second half of the extinction window.
#' @param minPE minimum preferred elevation of species, in m above sea level.
#' @param maxPE maximum preferred elevation of species, in m above sea level.
#' @param meanET mean elevation tolerance of species, in m.
#' @param sdET standard deviation of elevation tolerance of species, in m
#' @param meanLogPA mean base-10 logarithm of peak abundance of species.
#' @param sdLogPA standard deviation of the base-10 logarithm of peak abundance of  
#'   species.
#' @param maxPA maximum permitted peak abundance of species, as a percentage (e.g., set to  
#'   100 for the maximum peak abundance to be a 100% probability of collection)
#' @param minAff minimum affinity for channel facies (as opposed to floodplain facies).  
#'   Affinity can range from -1 to 1, with +1 for species preserved only in floodplain  
#'   facies, -1 for species preserved only in channel facies, and 0 for species preserved  
#'   equally in both facies.
#' @param maxAff maximum affinity for channel facies (as opposed to floodplain facies).  
#'   Affinity can range from -1 to 1, with +1 for species preserved only in floodplain  
#'   facies, -1 for species preserved only in channel facies, and 0 for species preserved  
#'   equally in both facies.
#' @param x,object an object of class `species`.
#' @param ... additional arguments to be passed.
#' 
#' @examples
#' spec <- species(timeStep=0.001, extRate=0.25, startingSpecies=1000, durationMy=5.0, 
#'   minPE=0, maxPE=200, meanET=10, sdET=2, meanLogPA=log(25), sdLogPA=log(5), maxPA=100)
#' print(spec)
#' plot(spec)
#' 
#' @rdname species
#' @export species
#'
#' @return species returns an object of class "species", which as print, summary, and plot methods.
#' 
#' A species object is a data frame with eight columns; each row corresponds to one species. `id` is a unique identifier corresponding to each species; this identifier starts at 1 and increases sequentially as each new species originates. `ancestor` corresponds to the id number of the species that gave rise to a species; its value is 0 for all of the seed species that are used to start the random-branching model. `origination` and `extinction` are the times of origination and extinction, in m.y. `PE` is the preferred environment (elevation, in m), that is, the elevation at which a species is most likely to occur. `ET` is the environmental tolerance (in m) reflecting the standard deviation of the species response curve; larger values correspond to species found over a greater range of elevations. `PA` is the peak abundance, corresponding to the probability of occurrence (expressed as a percentage; i.e., 100 indicates a 1.0 probability) of a species at its preferred elevation. `Aff` reflects the affinity a species has for being preserved in channel vs. floodplain facies.
#' 

species <- function(timeStep, extRate=0.25, startingSpecies=200, durationMy=10, peakExtFactor=1, peakExtDurationMy=1, peakExtTimeMy=2, minPE, maxPE, meanET, sdET, meanLogPA, sdLogPA, maxPA, minAff=-1, maxAff=1) {
	results <- randomBranchTimeVarying(timeStep, extRate, startingSpecies, durationMy, peakExtFactor, peakExtDurationMy, peakExtTimeMy)
	results <- appendEcologicalParameters(results, minPE, maxPE, meanET, sdET, meanLogPA, sdLogPA, maxPA, minAff, maxAff)
	class(results) <- "species"
	return(results)
}

#' @rdname species
#' @export

plot.species <- function(x, ...) {
	oldpar <- graphics::par(mfrow=c(3, 2), mar=c(5, 4, 0, 2) + 0.1)
	graphics::hist(x$origination, xlab='Time of origination (m.y.)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$extinction, xlab='Time of extinction (m.y.)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$PE, xlab='Preferred environment (m)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$ET, xlab='Environmental tolerance (m)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$PA, xlab='Peak Abundance (as %)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$Aff, xlab='Affinity', breaks=50, las=1, col='gray', main='')
	graphics::par(oldpar)
}

#' 
#' @rdname species
#' @export

print.species <- function(x, ...) {
	x
}

#' 
#' @rdname species
#' @export

summary.species <- function(object, ...) {
	cat("simulation duration:               ", max(object$extinction), "m.y.\n")
	cat("number of species:                 ", max(object$id), "\n")
	cat("range of preferred elevation (PE): ", min(object$PE), "-", max(object$PE), "m\n")
	cat("range of elevation tolerance (ET): ", min(object$ET), "-", max(object$ET), "m\n")
	cat("range of peak abundance (PA):      ", min(object$PA), "-", max(object$PA), "%\n")
	cat("range of affinity (Aff):           ", min(object$Aff), "-", max(object$Aff), "\n")
}
