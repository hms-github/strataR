#' @title Stratigraphic Columns
#'
#' @description Create Parameters for Carbonate Deposition.
#'
#' @details Create the parameters that describe how carbonate sediments will be deposited. These include a series of pairs of water depths and sedimentation rates, lag time, and initial water depth. Carbonate sedimentation rates are given by the pairs, with sedimentation rates at intermediate water depths given by linear interpolation. Lag time describes the length of time that carbonate sedimentation will be paused following the flooding of a subaerially exposed platform. Lag time has been found to be necessary to create the ubiquitous flooding surfaces of shallow-marine carbonates. Initial water depth specifies the water depth at the beginning of a simulation, allowing a simulation to begin at any water depth.
#'
#' @param waterDepths a vector containing the water depths (in m) at which sedimentation rate is specified. The first value must be 0.0, and subsequent values must increase monotonically.
#' @param rates a vector of the same length as waterDepths containing the corresponding sedimentation rates (in m / m.y.)
#' @param lagTime a numeric value giving the duration of time (in m.y.) that carbonate sedimentation will be paused following the flooding of a platform initially above sea level.
#' @param initialWaterDepth a numeric value giving the water depth at the beginning of the simulation of a stratigraphic column.
#'
#' @export
#' 
#' @return returns a list with lagTime (a vector with a single value), initialWaterDepth (a vector with a single value), and productionCurve (a data frame with pairs of values of water depth and sedimentation rate).
#'
#' @examples
#' depths <- c(0.0,  2.0,  5.0, 10.0, 40.0, 100.0)
#' rates  <- c(0.0, 40.0, 15.0, 10.0, 5.0, 2.0)
#' sedi <- carbSediment(depths, rates, lagTime=0.02, initialWaterDepth=2)
#' 

carbSediment <- function(waterDepths, rates, lagTime=0, initialWaterDepth=0) {
	if (min(waterDepths) > 0 | waterDepths[1] != 0) {
		warning("Minimum (first) water depth must be 0.0", call.=FALSE)
	}
	if (min(rates) > 0 | rates[1] != 0) {
		warning("Minimum (first) rate must be 0.0", call.=FALSE)
	}
	if (!all(waterDepths >= 0))  {
		warning("All water depths must be greater than or equal to zero", call.=FALSE)
	}
	if (! all(diff(waterDepths) > 0)) {
		warning("Water depths must be increasing order", call.=FALSE)
	}
	if (length(lagTime)>1) {
		warning("lagTime must be a single numeric value", call.=FALSE)
	}
	if (min(lagTime) < 0) {
		warning("lagTime must be non-negative", call.=FALSE)
	}
	if (length(initialWaterDepth)>1) {
		warning("initialWaterDepth must be a single numeric value", call.=FALSE)
	}
	productionCurve <- data.frame(waterDepth=waterDepths, rate=rates)
	
	carbSed <- list(lagTime=lagTime, initialWaterDepth=initialWaterDepth, productionCurve=productionCurve)
	
	class(carbSed) <- "carbSediment"
	return(carbSed)
}
