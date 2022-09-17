#' @title Plot First and Last Occurrences (FADs and LADs)
#'
#' @description Plot the number of first or last occurrences of specie through a stratigraphic column.
#'
#' @details The number of first or last occurrences is shown, binned over a specified stratigraphic interval. If a mass extinction is simulated, the window of the mass extinction and the horizon that corresponds to peak extinction rates are also indicated.
#'
#' @param occurrences an object of class [`occurrences`].
#' @param column an object of class [`column`].
#' @param type a string of either "fad" or "lad" specifying which is to be plotted.
#' @param sampleSpacing a numeric value indicating the binning interval (in meters) for  
#'   tallying first or last occurrences in the stratigraphic column.
#' @param xMax a numeric value for setting the largest number of first or last occurrences  
#'   to be shown. Use NA (the default) to not set a ceiling on this axis.
#' @param peakExtTimeMy a numeric value indicating when a mass extinction reached its peak  
#'   (in m.y.), if one was simulated. Use NA (the default) if no extinction was simulated.
#' @param extDurationMy a numeric value indicating when the duration of a mass extinction  
#'   (in m.y.), if one was simulated. Use NA (the default) if no extinction was simulated.
#' @param extinctionColor a color for displaying the window of extinction and the  
#'   stratigraphic position of peak extinction rates.
#'
#' @export
#' 
#' @examples
#' data(occu)
#' data(coluValley)
#' fadLadPlot(occurrences=occu, column=coluValley, type='fad', sampleSpacing=0.5)
#'

fadLadPlot <- function(occurrences, column, type=c('fad', 'lad'), sampleSpacing=0.5, xMax=NA, peakExtTimeMy=NA, extDurationMy=NA, extinctionColor='red') {
	fadsLads <- fadLad(occurrences)
	bins <- 0
	fadLadHist <- 0
	xlab <- ''
	
	if (type == 'fad') {
		upperLimit <- ceiling(max(fadsLads$fad))
		if (upperLimit %% sampleSpacing > 0) upperLimit <- upperLimit + sampleSpacing
		lowerLimit <- floor(min(fadsLads$fad))
		if (lowerLimit %% sampleSpacing > 0) lowerLimit <- 0 # base of section
		bins <- seq(lowerLimit, upperLimit, sampleSpacing)
		fadLadHist <- graphics::hist(fadsLads$fad, breaks=bins, plot=FALSE)
		xlab <- 'First occurrences'
	} else if (type == 'lad') {
		upperLimit <- ceiling(max(fadsLads$lad))
		if (upperLimit %% sampleSpacing > 0) upperLimit <- upperLimit + sampleSpacing
		lowerLimit <- floor(min(fadsLads$lad))
		if (lowerLimit %% sampleSpacing > 0) lowerLimit <- 0 # base of section
		bins <- seq(lowerLimit, upperLimit, sampleSpacing)
		fadLadHist <- graphics::hist(fadsLads$lad, breaks=bins, plot=FALSE)
		xlab <- 'Last occurrences'
	} else {
		stop('unrecognized plot type')
	}
	
	barHalfWidth <- sampleSpacing/2
	if (!is.na(xMax) &  xMax < max(fadLadHist$counts)) {
		fadLadHist$counts[fadLadHist$counts > xMax] <- xMax
	}
	xlim = range(fadLadHist$counts)
	if (!is.na(xMax)) {
		xlim = c(0, xMax)
	}
	ylim = c(min(column$stratPosition), max(column$stratPosition))
	plot(fadLadHist$counts, fadLadHist$mids, type='n', xlim=xlim, ylim=ylim, xlab=xlab, ylab='Stratigraphic Position (m)', las=1)
	
	if (massExtinctionOccurred(peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)) {
		addExtinctionBox(column=column, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy, leftEdge=min(fadLadHist$counts), rightEdge=max(fadLadHist$counts), extinctionColor=extinctionColor)
	}
	
	nonZero <- fadLadHist$counts > 0
	graphics::rect(0, fadLadHist$mids[nonZero]-barHalfWidth, fadLadHist$counts[nonZero], fadLadHist$mids[nonZero] + barHalfWidth, col='gray')
}
