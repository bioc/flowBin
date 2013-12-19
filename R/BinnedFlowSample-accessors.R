## =========================================================================
## "dumb" accessors: 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export
setGeneric("bin.labels", function(object) standardGeneric("bin.labels"))
#' @export
setGeneric("bin.labels<-", function(object,value) standardGeneric("bin.labels<-"))
#' @export
setMethod("bin.labels", signature="BinnedFlowSample", function(object) {object@bin.labels})
#' @export
setReplaceMethod("bin.labels", signature=c("BinnedFlowSample", "list"), function(object, value) {object@bin.labels <- value; object})

 

setGeneric("eventsInBins", function(object) standardGeneric("eventsInBins"))
# =========================================================================
#' Count number of events for each tube in each bin
#' @description Useful for QA of bin mapping
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @export
#' @aliases eventsInBins,BinnedFlowSample-method
setMethod(  "eventsInBins", 
            signature="BinnedFlowSample", 
            function(object) 
            {
               #First check that all bins were actually mapped.
               #Count number of unique bin labels in each tube, ignoring NA:
               num.bins <- sapply(bin.labels(object),function(x) length(unique(x[which(!is.na(x))])))
               if(length(unique(num.bins)) > 1)
                   stop('Mismatch in number of bins across tubes!')
                   

               countOneTube <- function(tube.labels)
               {
                    tube.unique <- sort(unique(tube.labels))
                    sapply(tube.unique, function(x) length(which(tube.labels==x)) )
               }
               bin.counts <- sapply(bin.labels(object),countOneTube)
               rownames(bin.counts) <- sort(unique(bin.labels(object)[[1]]))
               bin.counts
            })

