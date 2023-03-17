#' @title Basin Geometry
#'
#' @description Create a basin geometry object
#'
#' @details Creates a basin geometry object, which is needed to create subsidence, eustasy, sedimentation, and basin objects. A geometry object defines the spatial and temporal characteristics of a basin simulation. 
#' 
#' @param fallLineY vertical elevation of the fall line, in meters. Corresponds to the  
#'   edge of sedimentary basin, from which sediment is introduced.
#' @param shoreX distance of the initial shoreline from the fall line, in km.
#' @param deltaWidth width of the marine delta, in km.
#' @param deltaToeY initial water depth of the toe (base) of the marine delta, in m.
#' @param marginWidth width of the sedimentary basin, in km.
#' @param resolution a string (low, medium, or high) describing the spatial resolution of  
#'   the model. In most cases, medium is sufficient.
#' @param nonMarAlpha,marineAlpha a value from 0 to 1, describing the curvature of the  
#'   coastal plain and marine equilibrium profiles. Larger values create a more concave  
#'   profile
#' @param duration duration of the basin simulation, in m.y.
#' @param timeStep duration of one time step in the basin simulation, in m.y.
#' @param x,object an object of class `geometry`.
#' @param ... additional arguments to be passed.
#' 
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, 
#'   marginWidth=600, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' summary(geom)
#' 
#' @rdname geometry
#' @export geometry
#' 
#' @return geometry returns an object of class "geometry", which includes print and summary methods.
#'
#' A geometry object consists of a list of the twelve arguments used to create the geometry object. Vertical distances are in meters, horizontal distances are in kilometers, times are in millions of years. `nonMarAlpha` and `marineAlpha` are non-dimensional constants that describe the curvature of the nonmarine and marine equilibrium values. Larger values reflect a more strongly curved topographic profile.
#' 


geometry <- function(fallLineY, shoreX, deltaWidth, deltaToeY, marginWidth, resolution = c('medium', 'low', 'high'), nonMarAlpha, marineAlpha, duration, timeStep) {	
	resolution <- match.arg(resolution)
	
	if (missing(fallLineY)) {
		warning("Argument fallLineY must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(shoreX)) {
		warning("Argument shoreX must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(deltaWidth)) {
		warning("Argument deltaWidth must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(deltaToeY)) {
		warning("Argument deltaToeY must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (deltaToeY >= 0) {
		warning("Argument deltaToeY	must be negative (below sea level)", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(marginWidth)) {
		warning("Argument marginWidth must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(resolution)) {
		warning("Argument resolution must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(nonMarAlpha)) {
		warning("Argument nonMarAlpha must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(marineAlpha)) {
		warning("Argument marineAlpha must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(duration)) {
		warning("Argument duration must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (missing(timeStep)) {
		warning("Argument timeStep must be given a value", call.=FALSE, immediate.=TRUE)
	}
	
	if (resolution == 'high') {
		print("Using high resolution will greatly increase run times and file sizes, often at little gain over a low or medium resolution. Consider running a low-resolution or medium-resolution model first to see if that is sufficient.")
	}
	
	if (resolution == 'low') {
		print("Using low resolution can cause undesirable threshold effects. If the model plots reveal unexpected abrupt changes in sediment aggradation or partitioning, re-run the model at medium resolution.")
	}
	
	fallLineX <- 0     # by default, always at the left edge
	shoreY <- 0        # shore is defined to be at an elevation of zero
	deltaX <- switch(resolution, 'low'=1.0, 'medium'=0.1, 'high'=0.01)
	
	results <- list(fallLineX=fallLineX, fallLineY=fallLineY, shoreX=shoreX, shoreY=shoreY, deltaWidth=deltaWidth, marginWidth=marginWidth, deltaToeY=deltaToeY, deltaX=deltaX, nonMarAlpha=nonMarAlpha, marineAlpha=marineAlpha, duration=duration, timeStep=timeStep)
	class(results) <- "geometry"
	return(results)
}

#' 
#' @rdname geometry
#' @export

summary.geometry <- function(object, ...) {
	cat("Fall line elevation (fallLineY):              ", object$fallLineY, "m\n")
	cat("Shoreline placement (shoreX):                 ", object$shoreX, "km\n")
	cat("Delta width (deltaWidth):                     ", object$deltaWidth, "km\n")
	cat("Basin width (marginWidth):                    ", object$marginWidth, "km\n")
	cat("Depth of delta toe (deltaToeY):               ", object$deltaToeY, "m\n")
	cat("Spatial resolution (deltaX):                  ", object$deltaX, "m (set via resolution)\n")
	cat("Curvature of nonmarine profile (nonMarAlpha): ", object$nonMarAlpha, "(dimensionless)\n")
	cat("Curvature of marine profile (marineAlpha):    ", object$marineAlpha, "(dimensionless)\n")
	cat("Model duration (duration):                    ", object$duration, "m.y.\n")
	cat("Time step (timeStep):                         ", object$timeStep, "m.y.\n")
}
