## =========================================================================
## removeSparseBins -- remove bins with low numbers of events in one or more
##                     tubes. 
##               
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#Cutoff proportion is optional. If not specified, all populations with fewer than 2
#events are removed
removeSparseBins.def <- function(object, cutoff.prop=NULL)
{
    message('Filtering sparse bins.')
    #Before calling bin.counts, check if any bins failed to map to one of the tubes at all:
    all.bins <- Reduce(union, bin.labels(object))
    bins.correctly.mapped <- Reduce(intersect, bin.labels(object))
    to.remove <- setdiff(all.bins, bins.correctly.mapped)

    #First removal based on bin absence, prior to getting bin counts:
    removeLab <- function(tube.labels, lab.to.rem)
    {
        tube.labels[which(tube.labels==lab.to.rem)] <- NA
        tube.labels
    }
    for (lab.to.rem in to.remove)
    {
        bin.labels(object) <- lapply(bin.labels(object), removeLab, lab.to.rem)
    }

    bin.counts <- eventsInBins(object)
    #Decide which to remove based on <2:
    to.remove <- union(to.remove, names(which(apply(bin.counts, c(1), function(x){min(x) < 2 })))) #have to take names as index may be lost due to previous removal of bins
     
    #Decide which to remove based on cutoff:
    if(!is.null(cutoff.prop))
    {
        tube.totals <- sapply(tube.set(object), nrow)
        bin.props <- sapply(1:length(tube.set(object)), function(x){bin.counts[,x] / tube.totals[x]})
        to.remove.cutoff <- names(which(apply(bin.props, c(1), function(x){min(x) < cutoff.prop })))
        to.remove <- union(to.remove, to.remove.cutoff)
    }

    #Second removal based on bin counts:
    for (lab.to.rem in to.remove)
    {
        bin.labels(object) <- lapply(bin.labels(object), removeLab, lab.to.rem)
    }

   filtered.count <- eventsInBins(object) 
   mean.remaining <- mean(apply(filtered.count, c(2), sum))
   mean.before <- mean(sapply(tube.set(object), nrow))
   mean.removed <- mean.before - mean.remaining
   message(paste(length(to.remove), 'bins removed, containing a total of', format(mean.removed, digits=2), 'or', format(mean.removed/mean.before * 100, digits=2), '% of events (averaged across tubes).' ))

   object
}

#' Remove bins from a BinnedFlowSample with few events in them
#' 
#' @aliases removeSparseBins,BinnedFlowSample-method
#' 
#' @details This is important to do prior to calculating bin expression, as bins containing 2 or less events, for, example, cannot have their median computed.
#' 
#' @param object the BinnedFlowSample to act on
#' @param cutoff.prop=NULL the minimum proportion that a bin must contain to be kept. If NULL, only bins with no events in at least one tube will be removed.
#' 
#' @return a BinnedFlowSample with sparse bins removed
#' 
#'@export

setMethod('removeSparseBins', signature=signature("BinnedFlowSample"), definition=removeSparseBins.def )
