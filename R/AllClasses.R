# =========================================================================
#' A class similar to flowSet, but with extra information needed by flowBin
#' 
#' \describe{
#'	\item{\code{name}:}{\code{character} string - name of the object }
#'	   \item{\code{tube.set}:}{\code{list} of \code{flowFrame}s containing raw flow data. }
#'	   \item{\code{control.tubes}:}{Integer \code{vector} indicating which tubes in the list (if any) to use as negative controls. May be empty.}
#'	   \item{\code{bin.pars}:}{Integer \code{vector} indicating which parameters to use for binning. These must be in the same position in all tubes.}
#'	   \item{\code{measure.pars}:}{\code{list} of integer \code{vector}s indicating which parameters to use for measurement. These must be specified per tube.}
#'	   Note: all slots can be get and set using accessor methods, for example \code{bin.pars(myFlowSet) <- c(1,2,5)}
#' }
#
#' @export
#' @aliases name,FlowSample-method tube.set,FlowSample-method bin.pars,FlowSample-method 
#' bin.pars<-,FlowSample,vector-method  control.tubes,FlowSample-method 
#' control.tubes<-,FlowSample,vector-method measure.pars,FlowSample-method
#' measure.pars<-,FlowSample,list-method
#' name<-,FlowSample,character-method
#' tube.set<-,FlowSample,list-method
#' show,FlowSample-method
#' bin.pars bin.pars<-
#' control.tubes control.tubes<-
#' measure.pars measure.pars<-
#' name name<-
#' tube.set tube.set<-
#' FlowSample-class
#' @rdname FlowSample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("FlowSample", representation( name="character",
                                   tube.set="list",
                                   control.tubes="vector",
                                   bin.pars="vector",
                                   measure.pars="list"
                                   )#,
                               #  prototype = list(name=character(),
                               # tube.set=new('flowSet'),
                               #   control.tubes=vector(),
                               #   bin.pars=vector(),
                               #   measure.pars=vector()
                               # )
        )

# =========================================================================
#' A \code{FlowSample}, but with binning information for each tube
#' 
#' \describe{
#'	\item{\code{clust.labels}:}{\code{list} of cluster label vectors, one for each tube}
#'	   Note: all slots can be get and set using accessor methods, for example \code{bin.pars(myFlowSet) <- c(1,2,5)}
#' }
#' @export
#' @aliases bin.labels,BinnedFlowSample-method
#' bin.labels<-,BinnedFlowSample,list-method
#' bin.labels bin.labels<-
#' BinnedFlowSample-class
#' @rdname BinnedFlowSample
##
##  -  
## 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#`
setClass("BinnedFlowSample", representation(bin.labels="list"),
                             contains="FlowSample"
        )


