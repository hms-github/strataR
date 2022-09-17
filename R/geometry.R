#' Create a basin geometry object
#'
#' Creates a basin geometry object, which is needed to create subsidence, eustasy, sedimentation, and basin objects. A geometry object defines the spatial and temporal characteristics of a basin simulation.
#' 
#' @title geometry: The geometry function
#' @param fallLineY vertical elevation of the fall line, in meters. Corresponds to the edge of sedimentary basin, from which sediment is introduced.
#' @param shoreX distance of the initial shoreline from the fall line, in km.
#' @param deltaWidth width of the marine delta, in km.
#' @param deltaToeY initial water depth of the toe (base) of the marine delta, in m.
#' @param marginWidth width of the sedimentary basin, in km.
#' @param resolution a string (low, medium, or high) describing the spatial resolution of the model. In most cases, medium is sufficient.
#' @param nonMarAlpha a value from 0 to 1, describing the curvature of the coastal plain equilibrium profile.
#' @param marineAlphaa value from 0 to 1, describing the curvature of the marine equilibrium profile.
#' @param duration duration of the basin simulation, in m.y.
#' @param timeStep duration of one time step in the basin simulation, in m.y.
#' @examples
#' geom <- geometry(fallLineY=150, shoreX=200, deltaWidth=100, deltaToeY=-100, marginWidth=500, nonMarAlpha=0.5, marineAlpha=2.0, duration=3.0, timeStep=0.01)
#' 
#' @rdname geometry
#' @export geometry

geometry <- function(fallLineY, shoreX, deltaWidth, deltaToeY, marginWidth, resolution = c('medium', 'low', 'high'), nonMarAlpha, marineAlpha, duration, timeStep) {
	# horizontal distances (shoreX, deltaWidth, marginWidth) are in km
	# vertical distances (fallLineY, deltaToeY) are in m; positive is above sea level, negative is below
	# time values (duration, timeStep) are in m.y.
	# nonMarAlpha and marineAlpha are non-dimensional terms that describe the curvature of the alluvial plain and the marine profile
	
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

#' @return \code{NULL}
#' 
#' @rdname geometry
#' @export

plot.geometry <- function(x) {
	message(paste("A geometry cannot be plotted"))
}

#' @return \code{NULL}
#' 
#' @rdname geometry
#' @export

print.geometry <- function(x) {
	cat("fallLineX:   ", x$fallLineX, "km\n")
	cat("fallLineY:   ", x$fallLineY, "m\n")
	cat("shoreX:      ", x$shoreX, "km\n")
	cat("shoreY:      ", x$shoreY, "m\n")
	cat("deltaWidth:  ", x$deltaWidth, "km\n")
	cat("marginWidth: ", x$marginWidth, "km\n")
	cat("deltaToeY:   ", x$deltaToeY, "m\n")
	cat("deltaX:      ", x$deltaX, "m\n")
	cat("nonMarAlpha: ", x$nonMarAlpha, "(dimensionless)\n")
	cat("marineAlpha: ", x$marineAlpha, "(dimensionless)\n")
	cat("duration:    ", x$duration, "m.y.\n")
	cat("timeStep:    ", x$timeStep, "m.y.\n")
}

#' @return \code{NULL}
#' 
#' @rdname geometry
#' @export

summary.geometry <- function(x) {
	cat("Fall line elevation (fallLineY):              ", x$fallLineY, "m\n")
	cat("Shoreline placement (shoreX):                 ", x$shoreX, "km\n")
	cat("Delta width (deltaWidth):                     ", x$deltaWidth, "km\n")
	cat("Basin width (marginWidth):                    ", x$marginWidth, "km\n")
	cat("Depth of delta toe (deltaToeY):               ", x$deltaToeY, "m\n")
	cat("Spatial resolution (deltaX):                  ", x$deltaX, "m (set via resolution)\n")
	cat("Curvature of nonmarine profile (nonMarAlpha): ", x$nonMarAlpha, "(dimensionless)\n")
	cat("Curvature of marine profile (marineAlpha):    ", x$marineAlpha, "(dimensionless)\n")
	cat("Model duration (duration):                    ", x$duration, "m.y.\n")
	cat("Time step (timeStep):                         ", x$timeStep, "m.y.\n")
}
