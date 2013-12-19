checkQNorm.func <- function(object, normed.object, do.plot=TRUE)
{
    fp.model <- flowFPModel(flowSet(tube.set(object)), parameters=bin.pars(object), nRecursions=5)
    fp <- flowFP(flowSet(tube.set(object)), fp.model)
    
    par.mfrow <- par('mfrow')
    if(do.plot)
    {
    	par(mfrow=c(1,2))
      #pdf(paste(plot.dir,'/',name(object),'_prenorm.pdf', sep=''))
      plot(fp, type='qc') #make a switch plot=T for this?
      #dev.off()
    }
    fp.mat <- counts(fp, transformation = 'log2norm')
    sd.before <- apply(fp.mat, 1, sd) #return the QA data in the form of sds of each tube
    
    #normed.object <- quantileNormalise(object)
    fp.model <- flowFPModel(flowSet(tube.set(normed.object)), parameters=bin.pars(object), nRecursions=5)
    fp <- flowFP(flowSet(tube.set(normed.object)), fp.model)
    if(do.plot)
    {
      plot(fp, type='qc') #make a switch plot=T for this?\
    }
    par(mfrow =par.mfrow)
    fp.mat <- counts(fp, transformation = 'log2norm')
    sd.after <- apply(fp.mat, 1, sd)
    #sd.change <- sd.after - sd.before
    #before.mean <- mean(sd.before)
    #after.mean <- mean(sd.after)
    list(sd.before = sd.before, sd.after = sd.after)
}

#' @import flowFP
#' @name checkQNorm
#' @title Function to check the quantile normalisation of a FlowSample using flowFP
#' @details object and normed.object are compared using \code{flowFP} binning, to assess the deviation in bin counts between the two. 
#' @return \code{list} containing two matrices of standard deviations across bins (rows) vs tubes (columns) for before (\code{sd.before}) and after (\code{sd.after}).
#' @example inst/examples/quantileNormalise.R
#' @aliases checkQNorm, checkQNorm,FlowSample,FlowSample-method
#' @export
setMethod("checkQNorm", 
          signature=signature(object="FlowSample", normed.object="FlowSample"),
          definition=checkQNorm.func
          )