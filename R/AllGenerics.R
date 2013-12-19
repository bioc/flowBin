## AllGenerics.R -- all generic functions for package flowBin
## Also contains setMethod() definitions linking to functions
## Occasionally contains function itself when very short
## Accessors for most classes are kept in separate files.

## =========================================================================
## quantileNormalise -- 
## normalise to-bin parameters across all tubes
## @export
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


	#' @export
	setGeneric("quantileNormalise", function(object) standardGeneric("quantileNormalise"))

## =========================================================================
## checkQNorm -- return before and after SDs by flowFP for qnorm
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export checkQNorm
  setGeneric("checkQNorm", function(object, normed.object, do.plot=NULL) standardGeneric("checkQNorm"))

## =========================================================================
## mapBinsKNN -- use KNN to map given cluster labels for first tube
##                   in flowSample to other tubes. 
## 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


	#' @export
	setGeneric("mapBinsKNN", function(object, tube.1.labels, n.neighbours=NULL, snow.cluster=NULL, dequant=NULL) standardGeneric("mapBinsKNN"))


## =========================================================================
## kMeansBin -- use overfitted K means clustering to bin a sample
## KNN used to map clusters across tubes
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


	setGeneric("kMeansBin", function(object, n.bins=NULL, n.neighbours=NULL, snow.cluster=NULL, random.seed=NULL, dequantize=NULL) standardGeneric("kMeansBin"))


## =========================================================================
## kMeansBin -- use overfitted K means clustering to bin a sample
## KNN used to map clusters across tubes
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


	setGeneric("flowFPBin", function(object, n.bins=NULL, snow.cluster=NULL, dequantize=NULL) standardGeneric("flowFPBin"))



## =========================================================================
## getBinExpr -- extract expression levels of each bin, with options for
##               different methods, including using comparison to a negative
##               control
## 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


	setGeneric("getBinExpr", function(object, method=NULL, include.bin.medians=NULL, scale=NULL) standardGeneric("getBinExpr"))

## =========================================================================
## removeSparseBins -- remove bins with low numbers of events in one or more
##                     tubes. 
##               
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


	setGeneric("removeSparseBins", function(object, cutoff.prop=NULL) standardGeneric("removeSparseBins"))





## =========================================================================
## show Methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export
setMethod("show",
          signature=signature(object="FlowSample"),
          definition=function(object)
          {
              cat('A FlowSample containing', length(tube.set(object)), 'tubes\n')
              cat('\nControl tube(s): ', control.tubes(object), '\n')
              cat('Parameters to bin: ', bin.pars(object), '\n')
              cat('Parameters to measure: ')
              m.pars <- measure.pars(object)
              #If measure pars are the same for all tubes, just show the first:
              if(all(sapply(m.pars, function(x){length(x)})==length(m.pars[[1]])) && all(sapply(m.pars, function(x){x==m.pars[[1]]})))
              {
              	cat(str(m.pars[[1]], give.head=F))
              }
              else
              {
              	cat('\n')
              	for (i in 1:length(m.pars))
              	{
              		cat('Tube ',i, ': ')
              		cat(str(m.pars[[i]]))
              	}
              }
              
              	
          }
)
