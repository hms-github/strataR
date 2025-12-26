#' @title Calculated Thickness of Carbonate Deposition
#'
#' @description Calculate the thickness of carbonate deposition during a single time step when creating a stratigraphic column.
#'
#' @details Thickness of sediment is calculated.
#'
#' @param depth a number indicating the current water depth (in m).
#' @param timeStep a number indicating the duration (in m.y.) of each simulated time step.
#' @param carbSed an object of class `carbSediment`.
#'
#' @export
#' 
#' @return a numeric value of the thickness of carbonate sediment that was deposited.
#'
#' @examples
#' depths <- c(0.0,  2.0,  5.0, 10.0, 40.0, 100.0)
#' rates  <- c(0.0, 40.0, 15.0, 10.0, 5.0, 2.0)
#' carbSed <- carbSediment(depths, rates, lagTime=0.02, initialWaterDepth=2)
#' thickness <- carbThickness(depth=42, timeStep=0.01, carbSed=carbSed)
#' 

carbThickness <- function(depth, timeStep, carbSed) {
	carbProduction <- carbSed$productionCurve
	rate <- 0
# 	cat(paste("length of carb production:", nrow(carbProduction)))
	if (all(depth < carbProduction$waterDepth)) {          # above sea level
		rate <- carbProduction$rate[1]
		# if (is.null(rate)) {
# 			stop("Condition 1")
# 		}
	} else if (all(depth > carbProduction$waterDepth)) {   # very deep water
		rate <- carbProduction$rate[length(carbProduction$rate)]
		# if (is.null(rate)) {
# 			stop("Condition 2")
# 		}
	} else if (any(depth == carbProduction$waterDepth)) {  # depth is at a node
		rate <- carbProduction$rate[which(depth == carbProduction$waterDepth)]
		# if (is.null(rate)) {
# 			stop("Condition 3")
# 		}
	} else {                                               # depth is between two nodes
		deepNode    <- max(which(depth >= carbProduction$waterDepth))
		shallowNode <- min(which(depth <= carbProduction$waterDepth))
		deepDepth <- carbProduction$waterDepth[deepNode]
		shallowDepth <- carbProduction$waterDepth[shallowNode]
		deepRate <- carbProduction$rate[deepNode]
		shallowRate <- carbProduction$rate[shallowNode]
		position <- (depth - deepDepth) / (shallowDepth - deepDepth)
		rate <- deepRate + position * (shallowRate - deepRate)
		# if (is.null(rate)) {
# 			stop("Condition 4")
# 		}
	}
	thickness <- timeStep * rate
	thickness

	# prevent erosion below sea level or sedimentation above sea level
	if (abs(thickness) > abs(depth)) {
		thickness <- depth
	}
	return(thickness)
}
