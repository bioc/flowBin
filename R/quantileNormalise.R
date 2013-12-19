qnorm.func <- function(object)
{
    norm.sample <- object

    #Downsample to lowest cell count:
#     min.cells <- min(sapply(tube.set(object), function(x){nrow(exprs(x))}))
#     downSample <- function(this.frame, min.cells)
#     {
#         exprs(this.frame) <- exprs(this.frame)[sample(nrow(this.frame), min.cells),]
#         this.frame
#     }
    ##TODO! Put a QA check here that not losing too many cells
    ##TODO -- implement upSampling as well? Less lossy...
    #tube.set(norm.sample) <- lapply(tube.set(object), downSample, min.cells)

    normCol <- function(param, norm.sample)
    {
    	max.cells <- max(sapply(tube.set(norm.sample), function(x){nrow(exprs(x))}))
    	#Get the marker from each tube, filling in NAs to the length of the tube with most cells
        getCol <- function(the.tube)
        {
        	this.col <- exprs(the.tube)[,param]
        	this.col <- append(this.col, rep(NA, max.cells - length(this.col)))
        	
        	#Dequantise any margin events to stop them totally screwing quantile normalisation:
        	margin.values <- as.numeric(names(which(table(this.col) > nrow(the.tube)/10)))
        	for (margin.value in margin.values)
        	{
        		margin.events <- which(this.col==margin.value)
        		this.col[margin.events] <- this.col[margin.events] + runif(length(margin.events), 0, 0.0001)
        	}
        	
        	this.col
        }
    	
        param.mat <- sapply(tube.set(norm.sample), getCol)
        param.mat <- normalizeQuantiles(param.mat)
    	#Place the normed values back in:
        setCol <- function(tube.num)
        {
        	this.col <- param.mat[,tube.num]
        	this.col <- this.col[which(!is.na(this.col))]
        	exprs(tube.set(norm.sample)[[tube.num]])[,param] <- this.col
        	tube.set(norm.sample)[[tube.num]]
        }
        normed.tubes <- lapply(1:ncol(param.mat), setCol)
        tube.set(norm.sample) <- normed.tubes
        norm.sample
    }

    for(param in bin.pars(object))
        norm.sample <- normCol(param, norm.sample)

    norm.sample

}

## =========================================================================
#' quantileNormalise
#' normalise binning paramaters across all tubes of a flowSample
#' 
#' @description Since the binning parameters are the same across tubes, and samples
#' each tube is an aliquot from the same sample, these should have the same underlying
#' distribution. Hence, quantile normalisation can be used to force this to be so, 
#' removing technical variation.
#' @example inst/examples/quantileNormalise.R
#' 
#' @aliases quantileNormalise,FlowSample-method
#' 
#' @importFrom limma normalizeQuantiles
#' @export
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

setMethod("quantileNormalise", 
          signature=signature(object="FlowSample"),
          definition=qnorm.func
          )


