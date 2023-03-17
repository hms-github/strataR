#' @title Fossil Occurrences
#'
#' @description Generate occurrences of non-marine fossils in a stratigraphic column.
#'
#' @details Given a stratigraphic column and a set of species, `occurrences` will sample that column for the occurrence of species at a specified sample spacing. 
#'
#' @param column an object of class [`column`].
#' @param species an object of class [`species`].
#' @param sampleSpacing a number indicating the stratigraphic interval (in m) at which the  
#'   column will be sampled for fossil species.
#' @param x,object an object of class `occurrences`.
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
#' @return 'occurrences' returns an occurrences object, which includes print, summary, and plot methods.
#'
#' An occurrences object consists of a data frame with three columns; each row corresponds to one occurrence of a fossil species. `modelTime` is the age (in m.y.) of the occurrence, `stratPosition` is the position in meters above the base of column in which the species occurred, and `speciesId` is the id number of the species that occurred. This species id can be cross-referenced with the species object to find the times of origination and extinction, as well as the ecological and preservational characteristics of that species.
#'
#' @examples
#' data(coluValley)
#' data(spec)
#' occu <- occurrences(column=coluValley, species=spec, sampleSpacing=1.0)
#' summary(occu)
#' print(occu)
#' plot(occu, col=coluValley, species=spec)
#'
occurrences <- function(column, species, sampleSpacing=0.5) {
	numSpecies <- length(species$id)

	maxLength <- 100000
	storage <- vector(mode='numeric', length=maxLength)
	age       <- storage
	position  <- storage
	speciesId <- storage
	occElevation <- storage
	occFacies    <- storage
	
	collectionHorizons <- seq(0, max(column$stratPosition), sampleSpacing)
	occurrenceNum <- 1
	
	for (horizon in collectionHorizons) {
		horizonTime <- modelTimeForHorizon(column, horizon, oldestAge=FALSE)
		horizonElevation <- column$elevation[horizonTime == column$modelTime]
		horizonFacies <- column$facies[horizonTime == column$modelTime]
		for (i in 1:numSpecies) {
			if (horizonTime >= species$origination[i] && horizonTime <= species$extinction[i]) { # extant
				multiplier <- preservation(affinity=species$Aff[i], facies=horizonFacies)
				probCollection <- multiplier * probCollection(PE=species$PE[i], ET=species$ET[i], PA=species$PA[i], elevation=horizonElevation)						
				random <- stats::runif(1)
				if (random <= probCollection) {	                         # found
					age[occurrenceNum] <- horizonTime
					position[occurrenceNum] <- horizon
					speciesId[occurrenceNum] <- species$id[i]
					occElevation[occurrenceNum] <- horizonElevation
					occFacies[occurrenceNum] <- horizonFacies
					occurrenceNum <- occurrenceNum + 1
					if (occurrenceNum%%maxLength == 0) {
						age <- c(age, storage)
						position <- c(position, storage)
						speciesId <- c(speciesId, storage)
						occElevation<- c(horizonElevation, storage)
						occFacies <- c(horizonFacies, storage)
						message(paste("Resizing storage to hold up to", length(speciesId), "occurrences. If this happens multiple times, this may be a sign that the run needs to be terminated."))
					}
				}
			}
		}
	}	
	
	occurrences <- data.frame(modelTime=age, stratPosition=position, speciesId, elevation=occElevation, facies=occFacies)
	# trim any rows with no data
	occurrences <- occurrences[-which(occurrences$speciesId == 0), ]
	cat("Generated", nrow(occurrences), "occurrences\n")
	class(occurrences) <- "occurrences"
	return(occurrences)
}

#' 
#' @rdname occurrences
#' @export

plot.occurrences <- function(x, column, species, occurrenceColor='blue', rangeColor='blue', extantColor='gray', peakExtTimeMy=NA, extDurationMy=NA, extinctionColor='red', orderedByLad=TRUE, yAxisLabels=TRUE, ...) {
	# determine fads and lads
	fadsLads <- fadLad(x)
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
	
	plot(rows, fadsLads$fad, ylim=c(min(column$stratPosition), max(column$stratPosition)), xlab='Species', ylab=ylab, las=1, type='n', frame.plot=yAxisLabels, axes=yAxisLabels, ...)
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
		timeOrigination <- species$origination[theSpecies]
		timeExtinction  <- species$extinction[theSpecies]
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
#' @rdname occurrences
#' @export

print.occurrences <- function(x, ...) {
	x
}

#' 
#' @rdname occurrences
#' @export

summary.occurrences <- function(object, ...) {
	cat("oldest occurrence:   ", min(object$modelTime), "m.y.\n")
	cat("youngest occurrence: ", max(object$modelTime), "m.y.\n")
	cat("lowest occurrence:   ", min(object$stratPosition), "m\n")
	cat("highest occurrence:  ", max(object$stratPosition), "m\n")
	cat("occurring species:   ", length(unique(object$speciesId)), "\n")
}


