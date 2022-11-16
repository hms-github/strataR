#' @title Marine Fossil Occurrences
#'
#' @description Generate occurrences of marine fossils in a stratigraphic column.
#'
#' @details Given a stratigraphic column and a set of species, `marineOccurrences` will sample that column for the occurrence of species at a specified sample spacing.
#'
#' @param column an object of class [`column`].
#' @param marineSpecies an object of class [`marineSpecies`].
#' @param sampleSpacing a number indicating the stratigraphic interval (in m) at which the  
#'   column will be sampled for fossil species.
#' @param x,object an object of class `marineOccurrences`.
#' @param occurrenceColor,rangeColor,extantColor specify the color of the occurrence points, the preserved range of species, and the times in which the species was extant, including the time of origination and extinction.
#' @param peakExtTimeMy a numeric value (in m.y.) specifying the time of peak mass  
#'   extinction used; used in the occurrences plot.
#' @param extDurationMy a numeric value (in m.y.) specifying the duration of mass  
#'   extinction used; used in the occurrences plot.
#' @param extinctionColor a color specifying how a mass extinction interval should be  
#'   displayed on the occurrences plot.
#' @param orderedByLad a boolean specifying whether the occurrences should be ordered by  
#'   last occurrences (first occurrences if FALSE) on the occurrences plot.
#' @param yAxisLabels a boolean specifying whether y-axis labels should be displayed in  
#'   the occurrences plot.
#' @param ... additional arguments to be passed.
#'
#' @export
#' 
#' @return 'marineOccurrences' returns a marineOccurrences object, which includes print, summary, and plot methods.
#' 
#' An occurrences object consists of a list of three items, each row corresponding to one occurrence of a fossil species. `modelTime` is the age (in m.y.) of the occurrence, `stratPosition` is the position in meters above the base of column in which the species occurred, and `speciesId` is the id number of the species that occurred. This species id can be cross-referenced with the species object to find the times of origination and extinction, as well as the ecological and preservational characteristics of that species.
#'
#' @examples
#' data(marcolu)
#' data(marspec)
#' maroccu <- marineOccurrences(column=marcolu, marineSpecies=marspec, sampleSpacing=1.0)
#' summary(maroccu)
#' print(maroccu)
#' plot(maroccu, col=marcolu, marineSpecies=marspec)
#'
marineOccurrences <- function(column, marineSpecies, sampleSpacing=0.5) {
	# NOTE: because columns are in terms of waterDepth, anything that is below sea level will have a negative value, so must be converted to a positive value to convert to water depth	
	numSpecies <- length(marineSpecies$id)
	
	maxLength <- 100000
	storage <- vector(mode='numeric', length=maxLength)
	age       <- storage
	position  <- storage
	speciesId <- storage
	
	collectionHorizons <- seq(0, max(column$stratPosition), sampleSpacing)
	occurrenceNum <- 1
	
	for (horizon in collectionHorizons) {
		horizonTime <- modelTimeForHorizon(column, horizon, oldestAge=FALSE)
		horizonWaterDepth <- -column$elevation[horizonTime == column$modelTime]  # change sign to convert elevation to water depth
		for (i in 1:numSpecies) {
			if (horizonTime >= marineSpecies$origination[i] && horizonTime <= marineSpecies$extinction[i]) { # extant
				probCollection <- marineProbCollection(PD=marineSpecies$PD[i], DT=marineSpecies$DT[i], PA=marineSpecies$PA[i], waterDepth=horizonWaterDepth)						
				random <- stats::runif(1)
				if (random <= probCollection) {	                         # found
					age[occurrenceNum] <- horizonTime
					position[occurrenceNum] <- horizon
					speciesId[occurrenceNum] <- marineSpecies$id[i]
					occurrenceNum <- occurrenceNum + 1
					if (occurrenceNum%%maxLength == 0) {
						age <- c(age, storage)
						position <- c(position, storage)
						speciesId <- c(speciesId, storage)
						message(paste("Resizing storage to hold up to", length(speciesId), "occurrences. If this happens multiple times, this may be a sign that the run needs to be terminated."))
					}
				}
			}
		}
	}	
	
	occurrences <- as.data.frame(cbind(modelTime=age, stratPosition=position, speciesId))
	# trim any rows with no data
	occurrences <- occurrences[-which(occurrences$speciesId == 0), ]
	cat("Generated", nrow(occurrences), "occurrences")
	class(occurrences) <- "marineOccurrences"
	occurrences
}

#' 
#' @rdname marineOccurrences
#' @export

plot.marineOccurrences <- function(x, column, marineSpecies, occurrenceColor='blue', rangeColor='blue', extantColor='gray', peakExtTimeMy=NA, extDurationMy=NA, extinctionColor='red', orderedByLad=TRUE, yAxisLabels=TRUE, ...) {
	# determine fads and lads
	fadsLads <- marineFadLad(x)
	if (orderedByLad == TRUE) {
		fadsLads <- fadsLads[c(order(fadsLads$lad)),]
	} else {
		fadsLads <- fadsLads[c(order(fadsLads$fad)),]
	}
	rows <- seq(1:length(fadsLads$fad))
	
	ylab = 'Stratigraphic Position (m)'
	if (yAxisLabels == FALSE) {
		ylab = ''
	}
	
	plot(rows, fadsLads$fad, ylim=c(min(column$stratPosition), max(column$stratPosition)), xlab='Species', ylab=ylab, las=1, type='n', frame.plot=yAxisLabels, axes=yAxisLabels)
	if (yAxisLabels == FALSE) {
		graphics::axis(1)
	}
	
	# if there is a mass extinction, add a box for the interval and a line for the peak
	if (massExtinctionOccurred(peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy)) {
		addExtinctionBox(column=column, peakExtTimeMy=peakExtTimeMy, extDurationMy=extDurationMy, leftEdge=min(rows), rightEdge=max(rows), extinctionColor=extinctionColor)
	}
	
	# times of origination and extinction
	for (i in 1:length(fadsLads$occurringSpecies)) {
		theSpecies <- fadsLads$occurringSpecies[i]
		timeOrigination <- marineSpecies$origination[theSpecies]
		timeExtinction  <- marineSpecies$extinction[theSpecies]
		horizonOrigination <- stratPositionForAge(timeOrigination, column)
		horizonExtinction <- stratPositionForAge(timeExtinction, column)
		graphics::points(i, horizonOrigination, pch=3, cex=0.5, col=extantColor)
		graphics::points(i, horizonExtinction, pch=3, cex=0.5, col=extantColor)
		graphics::segments(i, horizonOrigination, i, horizonExtinction, lwd=0.5, lty='dotted', col=extantColor)
	}

	# add ranges
	graphics::segments(rows, fadsLads$fad, rows, fadsLads$lad, col=rangeColor, lwd=1)

	# add the occurrences to the plot
	for (i in 1:length(fadsLads$occurringSpecies)) {
		theSpecies <- fadsLads$occurringSpecies[i]
		horizons <- x$stratPosition[x$speciesId == theSpecies]
		speciesX <- rep(i, length(horizons))
		graphics::points(speciesX, horizons, col=occurrenceColor, pch=16, cex=0.7)
	}

	# singletons
	graphics::points(rows[fadsLads$numOccurrences == 1], fadsLads$fad[fadsLads$numOccurrences == 1], pch=16, cex=0.7, col='white')
	graphics::points(rows[fadsLads$numOccurrences == 1], fadsLads$fad[fadsLads$numOccurrences == 1], pch=1, cex=0.8, col=occurrenceColor)
}

#' 
#' @rdname marineOccurrences
#' @export

print.marineOccurrences <- function(x, ...) {
	cat("oldest occurrence:   ", min(x$modelTime), "m.y.\n")
	cat("youngest occurrence: ", max(x$modelTime), "m.y.\n")
	cat("lowest occurrence:   ", min(x$stratPosition), "m\n")
	cat("highest occurrence:  ", max(x$stratPosition), "m\n")
	cat("occurring species:   ", length(unique(x$speciesId)), "\n")
}

#' 
#' @rdname marineOccurrences
#' @export

summary.marineOccurrences <- function(object, ...) {
	cat("oldest occurrence:   ", min(object$modelTime), "m.y.\n")
	cat("youngest occurrence: ", max(object$modelTime), "m.y.\n")
	cat("lowest occurrence:   ", min(object$stratPosition), "m\n")
	cat("highest occurrence:  ", max(object$stratPosition), "m\n")
	cat("occurring species:   ", length(unique(object$speciesId)), "\n")
}


