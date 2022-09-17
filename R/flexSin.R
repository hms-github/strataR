#' @title Custom Sine Waves
#'
#' @description A method for customizing the shape of sine waves, for generating subsidence, eustatic sea level, and sediment flux histories.
#'
#' @details Change the amplitude and period of sine-wave curves. Specify one of four possible phases that describe the starting point of the curve, highPoint, lowPoint, or the falling or rising inflection points. Change how flattened the curves are, ranging from classical sine waves to square waves. Also change the symmetry of the curves, with symmetrical rises and falls, or with one side being steeper than the other.
#'
#' @param x a numeric vector of time points, in millions of years.
#' @param period period of cosine function, in millions of years.
#' @param amplitude amplitude of the cosine function, in meters. Note that this is half  
#'   the peak-to-peak height of the curve.
#' @param symmetry a number from 0 to 1 describing the symmetry of the curve. 0.5 is a  
#'   perfectly symmetrical wave. Values approaching zero have increasingly faster and  
#'   shorter falls and slower longer rises. Values approaching one have the opposite:  
#'   increasingly slow and longer falls and faster and shorter rises.
#' @param shape a positive number reflecting how squared the curve is. 0 corresponds to a  
#'   classical cosine curve, and increasingly positive values produce a wave that is  
#'   increasingly square-shaped.
#' @param phase a string of one of four values describing whether the curve starts at the  
#'   high point ("highPoint"), low point ("lowPoint"), at the inflection point of the  
#'   rising limb ("rising"), or at the inflection point of the falling limb ("falling").
#'
#' @export
#' 
#' @return a sine curve of a specified shape
#'
#' @examples
#' timePoints <- seq(0.0, 10.0, 0.001) # 10 million years in 1000-year steps
#' sineWave <- flexSin(timePoints, period=3, amplitude=10, symmetry=0.2, shape=4, 
#'   phase="rising")
#' plot(timePoints, sineWave, type="l")
#'

flexSin <- function(x, period=1, amplitude=1, symmetry=0.5, shape=0, phase=c("falling", "rising", "highPoint", "lowPoint")) {
	midpoint <- 0.5
	position <- (x/period)%%1
	position <- phaseScaling(position, phase=phase, symmetry=symmetry)
	position <- sapply(position, symmetryScaling, symmetry=symmetry)
	y <- amplitude * sqrt((1+shape^2) / (1+shape^2*cos(position * 2 * pi)^2)) * cos(position * 2 * pi)
	y <- y - y[1]  # shift curve so that initial y-value is always zero (starts at origin)
	y
	return(y)
}