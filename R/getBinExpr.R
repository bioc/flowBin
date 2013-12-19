## =========================================================================
## Various distance metric functions to use with getBinExpr: 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

medianFI <- function(test, control)
{
    median(test)
}

#Distance function: difference of MFIs
medianFIDist <- function(test, control)
{
    if(is.null(control))
      stop('NULL control frame -- no control tubes specified?')
    #Note: assumes log transform has been applied, so relinearises before subtracting
    res <- median(10^test) - median(10^control)
    if(res < 1) res <- 1
    log(res,10)
}

proportionPositive <- function(test, control)
{
	if(is.null(control))
		stop('NULL control frame -- no control tubes specified?')
	
	threshold <- quantile(control,0.98)
	res <- length(which(test > threshold)) / length(test)
	res
}

## =========================================================================
# getBinExpr main function definition 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

.getBinExpr.def <- function(object, method='medianFI', include.bin.medians=T, scale=T)
{
    message(paste('Calculating', method ,'for all populations.'));
    bins <- unique(bin.labels(object)[[1]]) 
    bins <- bins[!is.na(bins)]
    distance.func <- switch(method, 
                            medianFI=medianFI,
                            medianFIDist=medianFIDist,
    						propPos=proportionPositive
                           )
    #Function to get a frame containing only events in the given bin
    getBinFrame <- function(the.frame, bin.labels, bin.to.get)
    {
        if(length(which(bin.labels==bin.to.get)) < 2)
            stop(paste("Bin ", bin.to.get, "has less than 2 events in at least one tube! Consider using removeSparseBins."))
        bin.frame <- the.frame
        exprs(bin.frame) <- exprs(the.frame)[which(bin.labels==bin.to.get),]
        bin.frame
    }

    #Function to scale extracted data using range as set in FCS:
    scaleData <- function(par.col, data, object )
    {
        par.min <- pData(parameters(object@tube.set[[1]]))[par.col,'minRange']
        par.max <- pData(parameters(object@tube.set[[1]]))[par.col,'maxRange']
        (data[,as.character(par.col)] - par.min) / (par.max - par.min)
    }

    #Function to get expr of a single population in terms of all measured parameters across all tubes:
    getOnePopExpr <- function(object, pop.label, distance.func)
    {
        #Add code to get ordering of tubes if named:
        #if(!is.numeric(control.tubes(object)))
        #    this.control.tubes <- which(names())
        
        test.tubes <- setdiff(1:length(tube.set(object)), control.tubes(object))      
        extractOneTube <- function(tube.no, object)
        {
            
            #Use the closest control before the current tube (if multiple controls given) 
            if(length(control.tubes(object))>0)
            {
              which.control <- control.tubes(object)[max(which(tube.no - control.tubes(object) > 0))]
              control.frame <- getBinFrame(tube.set(object)[[which.control]], bin.labels(object)[[which.control]], pop.label )
            }
            else
            {
              #If no control tubes, return a NULL control frame
              control.frame <- NULL
            }
            
            test.frame <- getBinFrame(tube.set(object)[[tube.no]], bin.labels(object)[[tube.no]], pop.label )

            test.params <- measure.pars(object)[[tube.no]]

            tube.expr <- sapply(test.params, function(test.param){distance.func(exprs(test.frame)[,test.param], exprs(control.frame)[,test.param])})
            if(scale==T&&method!='propPos')
            {
                tmp.expr <- matrix(tube.expr, ncol=length(tube.expr))
                colnames(tmp.expr) <- test.params
                scaled.expr <- sapply(test.params, scaleData, tmp.expr, object)
                tube.expr <- scaled.expr 
                tube.expr[which(tube.expr < 0)] <- 0
            }
            names(tube.expr) <- pData(parameters(test.frame))$desc[test.params]
            tube.expr
        }
        Reduce(`c`, lapply(test.tubes, extractOneTube, object))
    }
   
    getOneBinMeds <- function(object, pop.label)
    {
      if(length(control.tubes(object))>0)
      {
        which.tube <- control.tubes(object)[1]
      }else
      {
        which.tube <- setdiff(1:length(tube.set(object)), control.tubes(object))[1]
      }
      the.frame <- getBinFrame(tube.set(object)[[which.tube]], bin.labels(object)[[which.tube]], pop.label )
      bin.meds <- sapply(bin.pars(object),function(x) median(exprs(the.frame)[,x]) )
                
      names(bin.meds) <- pData(parameters(the.frame))$desc[bin.pars(object)]

      bin.meds

    }

    #Now apply over all bins to get expression data:
    all.expr <- sapply(bins, function(x,y){
      #message(paste('Calculating', method ,'for bin',x)); 
      getOnePopExpr(y,x, distance.func)
      }, object)
    all.expr <- t(all.expr)

    #Add medians from binned channels if requested:
    if(include.bin.medians==T)
    {
    	#NB: if everything breaks, it's from this.
        all.meds <- matrix(sapply(bins, function(x,y){getOneBinMeds(y,x)}, object), nrow=length(bin.pars(object)))
        all.meds <- t(all.meds)
        #Scale bin medians if requested:
        if(scale==T)
        {
            tmp.meds <- all.meds
            colnames(tmp.meds) <- bin.pars(object)
            scaled.meds <- sapply(bin.pars(object), scaleData, tmp.meds, object)
            colnames(scaled.meds) <- colnames(all.meds)      
            all.meds <- scaled.meds
        }
        all.expr <- cbind(all.meds, all.expr)
    }

    rownames(all.expr) <- paste(name(object), bins, sep=' pop ')
    all.expr
    
    #Also will need code to deal with weird cases like duplicate parameters
}

#' getBinExpr main function definition
#' @name getBinExpr
#' @title Calculate the expression of each bin in a BinnedFlowSample in terms of the measurement markers
#' 
#' @aliases getBinExpr getBinExpr,BinnedFlowSample-method
#' @param method Method to use to compute expression, passed as a string. Defaults to \code{medianFI}, which takes the simple median of each bin,
#' and does not require control tubes. Other options available are \code{medianFIDist}, which uses medians with the 
#' median of the negative control subtracted out, and \code{propPos} which sets a threshold at the 98th percentile of the
#' negative control and determines what proportion of cells lie above that.
#' @param include.bin.medians logical, specifies whether to compute the medians of each bin in terms of the binning
#' markers and include them in the result or not. Defaults to \code{T}.
#' @param scale logical specifying whether to scale the results to the interval (0,1). If \code{T} (default), 
#' then all medians will be divided by the range for that marker as specified in the flowFrame.
#' 
#' @example inst/examples/getBinExpr.R
#' 
#' @return A numeric matrix containing expression values, with bins as rows and markers as columns
#' 
#' @export
setMethod("getBinExpr", signature=signature("BinnedFlowSample"), definition=.getBinExpr.def)

