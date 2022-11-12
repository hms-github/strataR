#' @title Marine Species Simulation
#'
#' @description Create a marine species object.
#'
#' @details Creates a marine species object, which records the ancestor of a species, its times of origination and extinction, and its ecological characteristics: preferred water depth, water depth tolerance, and peak abundance. Evolution of species is treated as a random-branching process, which can be optionally be decreased during a mass extinction. A marine species object is designed for marine species; for nonmarine species, see species.
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
#' @param minPD minimum preferred water depth of species, in m below sea level.
#' @param maxPD maximum preferred water depth of species, in m below sea level.
#' @param meanDT mean water depth tolerance of species, in m.
#' @param sdDT standard deviation of water depth tolerance of species, in m.
#' @param meanLogPA mean base-10 logarithm of peak abundance of species.
#' @param sdLogPA standard deviation of the base-10 logarithm of peak abundance of  
#'   species.
#' @param maxPA maximum permitted peak abundance of species, as a percentage (e.g., set to  
#'   100 for the maximum peak abundance to be a 100% probability of collection).
#' @param x,object an object of class `marineSpecies`.
#' @param ... additional arguments to be passed to the plot.
#' 
#' @examples
#' spec <- marineSpecies(timeStep=0.001, extRate=0.25, startingSpecies=1000, 
#'   durationMy=5.0, minPD=0, maxPD=200, meanDT=10, sdDT=2, meanLogPA=log(25), 
#'   sdLogPA=log(5), maxPA=100)
#' summary(spec)
#' print(spec)
#' plot(spec)
#' 
#' @rdname marineSpecies
#' @export marineSpecies
#' 
#' @return marineSpecies returns an object of class "marineSpecies", which as print, summary, and plot methods.
#'
#' A marine species object consists of a data frame with seven columns; each row corresponds to one species. `id` is a unique identifier corresponding to each species; this identifier starts at 1 and increases sequentially as each new species originates. `ancestor` corresponds to the id number of the species that gave rise to a species; its value is 0 for all of the seed species that are used to start the random-branching model. `origination` and `extinction` are the times of origination and extinction, in m.y. `PD` is the preferred depth (water depth, in m), that is, the water depth at which a species is most likely to occur. `DT` is the depth tolerance (in m) reflecting the standard deviation of the species response curve; larger values correspond to species found over a greater range of water depths. `PA` is the peak abundance, corresponding to the probability of occurrence (expressed as a percentage; i.e., 100 indicates a 1.0 probability) of a species at its preferred water depth.
#' 

marineSpecies <- function(timeStep, extRate=0.25, startingSpecies=200, durationMy=10, peakExtFactor=1, peakExtDurationMy=1, peakExtTimeMy=2, minPD, maxPD, meanDT, sdDT, meanLogPA, sdLogPA, maxPA) {
	results <- randomBranchTimeVarying(timeStep, extRate, startingSpecies, durationMy, peakExtFactor, peakExtDurationMy, peakExtTimeMy)
	results <- marineAppendEcologicalParameters(results, minPD, maxPD, meanDT, sdDT, meanLogPA, sdLogPA, maxPA)
	class(results) <- "marineSpecies"
	return(results)
}

#' @rdname marineSpecies
#' @export

plot.marineSpecies <- function(x, ...) {
	oldpar <- graphics::par(mfrow=c(3, 2), mar=c(5, 4, 0, 2) + 0.1)
	graphics::hist(x$origination, xlab='Time of origination (m.y.)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$extinction, xlab='Time of extinction (m.y.)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$PD, xlab='Preferred depth (m)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$DT, xlab='Depth tolerance (m)', breaks=50, las=1, col='gray', main='')
	graphics::hist(x$PA, xlab='Peak Abundance (as %)', breaks=50, las=1, col='gray', main='')
	plot(1, 1, type='n', axes=FALSE, xlab="", ylab="", ...)  # blank for remaining cell
	graphics::par(oldpar)
}

#' 
#' @rdname marineSpecies
#' @export

print.marineSpecies <- function(x, ...) {
	cat("simulation duration:           ", max(x$extinction), "m.y.\n")
	cat("number of species:             ", max(x$id), "\n")
	cat("range of preferred depth (PD): ", min(x$PD), "-", max(x$PD), "m\n")
	cat("range of depth tolerance (DT): ", min(x$DT), "-", max(x$DT), "m\n")
	cat("range of peak abundance (PA):  ", min(x$PA), "-", max(x$PA), "%\n")
}

#' 
#' @rdname marineSpecies
#' @export

summary.marineSpecies <- function(object, ...) {
	cat("simulation duration:           ", max(object$extinction), "m.y.\n")
	cat("number of species:             ", max(object$id), "\n")
	cat("range of preferred depth (PD): ", min(object$PD), "-", max(object$PD), "m\n")
	cat("range of depth tolerance (DT): ", min(object$DT), "-", max(object$DT), "m\n")
	cat("range of peak abundance (PA):  ", min(object$PA), "-", max(object$PA), "%\n")
}
