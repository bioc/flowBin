#' @docType data
#' @title Multitube AML sample as example data for flowBin
#' @keywords datasets
#' @name amlsample
#' @usage data(amlsample)
#' @aliases aml.sample
#' @format a \code{flowSample} containing 7 tubes with 3 common parameters and 4 measure parameters per tube.
#' @import methods
#' @source FlowRepository.org accession FR-FCM-ZZYA

NULL



## =========================================================================
#' function to run the entire flowBin pipeline
#' @aliases flowbin-package
#' @description Takes a list of flowFrames representing tubes from a single flow cytometry sample, and combines them using binning of events in terms of common markers across tubes.
#' @param tube.list a list of flowFrames, one for each tube to combine
#' @param bin.pars a numerical vector indicating which flow parameters in the 
#' each flowFrame to use for combining tubes. These should be the same markers
#' assayed across every tube.
#' @param control.tubes a vector indicating which tubes in tube.list to use for negative controls. May be empty.
#' @param measure.pars a list of which parameters to measure expression for, 
#' with one vector for each tube. If left NULL, this defaults to all parameters
#' other than those specified as \code{bin.pars}
#' @param sample.name name of this flowSample, for convenience (defaults to 'Unnamed Flow Expr Set')
#' @param bin.method The method to use for creating bins. The two options are "kmeans" for k-means 
#' clustering and nearest-neighbour mapping of bins. or "flowFP" for flowFP binning and direct
#' mapping of bin boundaries across tubes.
#' @param expr.method The method to use to compute bin expression across tubes. This defaults to 
#' MFI of the cells belonging to that bin in each tube. Other options are 
#' @param sparse.bin.thresh Bins which contain fewer than this proportion of total events in any tube will be excluded as outliers. Defaults to 0.001
#' @param dequantize If TRUE, adds a small (region of 1e-8) value to flow data to help break ties when binning.
#' @param snow.cluster A cluster created using the \code{snow} package, which flowBin will use to speed up computation. If NULL, flowBin will execute in serial mode.
#' @param n.bins Number of bins to use. Note that this must be a power of 2 if flowFP is selected as binning method.
#' @param scale.expr If TRUE, the resulting expression values will be scaled to (0,1) using the ranges specified in the flowFrames in \code{tube.list}.
#' @param do.qnorm If TRUE, the binning markers will be quantile normalized prior to binning.
#' @param return.bins If TRUE, return a \code{BinnedFlowExprSet} containing the bins themselves as well as the expression for each bin.
#' @return A matrix containing expression values for each bin in terms of each marker across all tubes. If \code{return.bins} is set TRUE, then a list containing a \code{BinnedFlowExprSet} followed by the expression matrix is returned.
#' @export
#' @importFrom snow clusterEvalQ
#' @example inst/examples/flowBin.R
#' 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

flowBin <- function(tube.list, 
                    bin.pars, 
                    control.tubes=vector(), 
                    measure.pars=NULL, 
                    sample.name='Unnamed Flow Expr Set', 
                    bin.method='kmeans',
                    expr.method='medianFI', 
                    sparse.bin.thresh=0.001, 
                    dequantize=T, 
                    snow.cluster=NULL, 
                    n.bins=128, 
                    scale.expr=F, 
                    do.qnorm=T,
					return.bins=F
                    )
{
  ####Argument processing ####
  
  #If no measure.pars given, use all pars except bin.pars
  if(is.null(measure.pars))
    measure.pars <- lapply(tube.list, function(x){setdiff(1:ncol(x), bin.pars)})
     
  #####Make flowSet and run flowBin ####
  this.sample <- new('FlowSample', name=sample.name, tube.set=tube.list, control.tubes=control.tubes, bin.pars=bin.pars, measure.pars=measure.pars)
  message(paste('Applying flowBin to', sample.name))
  
  if(do.qnorm)
  {
  	message('Quantile normalising binning parameters across tubes')
    this.sample <- quantileNormalise(this.sample)
  }
  
  #TODO: if flowFP....
  if(!bin.method %in% c('kmeans', 'flowFP'))
  	stop('Bin method not recognised. It should be kmeans or flowFP.')
  message(paste('Binning using', bin.method))
  if(bin.method=='kmeans')
  {
  	binned.sample <- kMeansBin(this.sample, dequantize=dequantize, n.bins=n.bins)
  	binned.sample <- removeSparseBins(binned.sample, sparse.bin.thresh) #TODO: empirically determine best default
  }
  
  if(bin.method=='flowFP')
  {
  	binned.sample <- flowFPBin(this.sample, n.bins=n.bins, dequantize=dequantize)
  	binned.sample <- removeSparseBins(binned.sample, sparse.bin.thresh) #TODO: empirically determine best default
  }
  sample.expr <- getBinExpr(binned.sample, method=expr.method, scale=scale.expr)
  if(return.bins)
  	return(list(binned.sample=binned.sample, sample.expr=sample.expr))
  
  return(sample.expr)
  
}