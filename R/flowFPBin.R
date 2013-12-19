flowFPBin.def <- function(object, n.bins=128, snow.cluster=NULL, dequantize=T)
{
  	n.div <- floor(log(n.bins,2))
  	
    #First, turn frame list into flowSet
  	tryCatch(the.fs <- flowSet(object@tube.set), error = function(e) stop('Cannot convert this data into a uniform flowSet for flowFP. Try another binning method.'))
  	#If error, can't operate on it, sorry.
  		
  	#Next, train and microbin
  	the.model <- flowFP::flowFPModel(the.fs, parameters=bin.pars(object), nRecursions=n.div)
      
    bin.labels <- lapply(object@tube.set, function(the.frame){flowFP::flowFP(the.frame, the.model)@tags[[1]]})
  	
  	new('BinnedFlowSample', bin.labels=bin.labels, 
  		name=name(object), 
  		tube.set=tube.set(object), 
  		control.tubes=control.tubes(object), 
  		bin.pars=bin.pars(object), 
  		measure.pars=measure.pars(object))
}

## =========================================================================
#' Bin sample using flowFP binning
#' @aliases flowFPBin,FlowSample-method
#' 
#' @param object flowSample to bin
#' @param n.bins=128 number of bins to use. This should be a power of 2, and will be rounded down to the nearest power of 2 if not.
#' @param snow.cluster=NULL Optional snow cluster to use for parallel execution.
#' @param dequantize=T If TRUE, adds a small (region of 1e-8) value to flow data to help break ties when binning.
#' 
#' @example inst/examples/flowFPBin.R
#' 
#' @return a \code{BinnedFlowSample}
#' @export
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

setMethod("flowFPBin", signature=signature("FlowSample"), definition=flowFPBin.def)
