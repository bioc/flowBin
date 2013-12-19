mapBinsKNN.def <- function(object, tube.1.labels, n.neighbours=1, snow.cluster=NULL, dequant=T)
{


    mapClustersKNN <- function(frame.dest, frame.source, cl, n.neighbours, bin.p, dequant)
    {
        source.exprs <- matrix(exprs(frame.source)[, bin.p], ncol=length(bin.p))
        dest.exprs <- matrix(exprs(frame.dest)[, bin.p], ncol=length(bin.p))
        #Ramp dequantisation method borrowed from flowFP:
        dequant.alpha <- 1e-8
        if(dequant)
        {
	    ramp.source <- seq(dequant.alpha, dequant.alpha * nrow(source.exprs), by=dequant.alpha)
            source.exprs <- source.exprs + ramp.source
	    ramp.dest <- seq(dequant.alpha, dequant.alpha * nrow(dest.exprs), by=dequant.alpha)
            dest.exprs <- dest.exprs + ramp.dest
        }
        as.integer(knn(source.exprs, dest.exprs, cl, n.neighbours))
    }

    #Util function to turn a flowSet back into a list to use SNOW apply:
    fsToList <- function(the.set)
    {
        ret.list <- list()
        for(i in 1:length(the.set))
            ret.list[[i]] <- the.set[[i]]
        ret.list
    }

    #If snow cluster provided, do KNN in parallel, otherwise do it serially:
    if(!is.null(snow.cluster))
    {     
        dump <- capture.output(clusterEvalQ(snow.cluster, library(flowBin)))
        #dump <- capture.output(clusterEvalQ(snow.cluster, require(class)))

        rest.labels <- parLapply(snow.cluster, fsToList(tube.set(object)[-1]), mapClustersKNN, tube.set(object)[[1]], tube.1.labels, n.neighbours, bin.pars(object), dequant)

    } else
    {
        #Note: serial not yet tested! Try out with smaller data set:
        rest.labels <- lapply(tube.set(object)[-1], mapClustersKNN, tube.set(object)[[1]], tube.1.labels, n.neighbours, bin.pars(object), dequant)
    }

    clust.labels <- append(list(tube.1.labels), rest.labels)
    new('BinnedFlowSample', bin.labels=clust.labels, 
                            name=name(object), 
                            tube.set=tube.set(object), 
                            control.tubes=control.tubes(object), 
                            bin.pars=bin.pars(object), 
                            measure.pars=measure.pars(object))
}

#' Internal function to map bins by KNN
#' @aliases mapBinsKNN,FlowSample-method
#' 
#' @details Takes a FlowSample and labels for the events in tube 1, and maps these to all other tubes.
#' 
#' @param object flowSample to map the bins of
#' @param tube.1.labels integer vector of bin labels for the events in tube 1
#' @param n.neighbours=1 number of neighbours to use for KNN mapping of bins from clustered tube
#' @param snow.cluster=NULL Optional snow cluster to use for parallel execution.
#' @param dequant=T If TRUE, adds a small (region of 1e-8) value to flow data to help break ties when binning.
#' 
#' @return a \code{BinnedFlowSample}
#' 
#' @example inst/examples/mapBinsKNN.R
#' 
#' @importFrom class knn
#' @importFrom BiocGenerics counts
#' @export
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## =========================================================================

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


setMethod("mapBinsKNN", signature=signature("FlowSample"), definition=mapBinsKNN.def)

