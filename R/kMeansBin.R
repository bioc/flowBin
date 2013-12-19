kMeansBin.def <- function(object, n.bins=100, n.neighbours=1, snow.cluster=NULL, random.seed=101, dequantize=T)
{
    set.seed(random.seed)
    #Get to-bin matrix
    #If only one bin par, R returns a vector instead of a matrix, so turn it back to a matrix:
    
    to.bin <- matrix(exprs(tube.set(object)[[1]])[, bin.pars(object)], ncol=length(bin.pars(object)))
    
    
    
    
    #Scale: -- get min, max for each parameter, scale as (data - min) / (max - min)
    this.pd <- pData(parameters(tube.set(object)[[1]]))
    rownames(this.pd) <- this.pd$name
    for( bin.par in 1:length(bin.pars(object)))
    {
    	bin.col <- bin.pars(object)[bin.par]
    	to.bin[ ,bin.par] <- (to.bin[ ,bin.par] - this.pd[bin.col,'minRange']) / (this.pd[bin.col,'maxRange'] - this.pd[bin.col, 'minRange'])
    }

    #Cluster:
    tube1.labels <- kmeans(to.bin, n.bins, iter.max=100)$cluster

    mapBinsKNN(object, tube1.labels, n.neighbours, snow.cluster, dequantize)
}

#' Bin sample using K-means binning
#' @aliases kMeansBin,FlowSample-method
#' 
#' @details Runs K-means clustering on the binning markers in the first tube of the data set. 
#' These clusters are then mapped to the other tubes using K-nearest neighbours.
#' 
#' @param object flowSample to bin
#' @param n.bins=128 number of bins to use. This should be a power of 2, and will be rounded down to the nearest power of 2 if not.
#' @param n.neighbours=1 number of neighbours to use for KNN mapping of bins from clustered tube
#' @param snow.cluster=NULL Optional snow cluster to use for parallel execution.
#' @param random.seed=101 Random seed to set to make K-means clustering deterministic.
#' @param dequantize=T If TRUE, adds a small (region of 1e-8) value to flow data to help break ties when binning.
#' 
#' @example inst/examples/kMeansBin.R
#' 
#' @return a \code{BinnedFlowSample}
#' 
#' @importFrom class knn
#' @export
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

setMethod("kMeansBin", signature=signature("FlowSample"), definition=kMeansBin.def)
