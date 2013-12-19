\name{FlowSample}
\alias{bin.pars}
\alias{bin.pars<-}
\alias{bin.pars,FlowSample-method}
\alias{bin.pars<-,FlowSample,vector-method}
\alias{control.tubes}
\alias{control.tubes<-}
\alias{control.tubes,FlowSample-method}
\alias{control.tubes<-,FlowSample,vector-method}
\alias{FlowSample}
\alias{FlowSample-class}
\alias{measure.pars}
\alias{measure.pars<-}
\alias{measure.pars<-,FlowSample,list-method}
\alias{measure.pars,FlowSample-method}
\alias{name}
\alias{name<-}
\alias{name<-,FlowSample,character-method}
\alias{name,FlowSample-method}
\alias{show,FlowSample-method}
\alias{tube.set}
\alias{tube.set<-}
\alias{tube.set<-,FlowSample,list-method}
\alias{tube.set,FlowSample-method}
\title{A class similar to flowSet, but with extra information needed by flowBin}
\description{
  \describe{ \item{\code{name}:}{\code{character} string -
  name of the object } \item{\code{tube.set}:}{\code{list}
  of \code{flowFrame}s containing raw flow data. }
  \item{\code{control.tubes}:}{Integer \code{vector}
  indicating which tubes in the list (if any) to use as
  negative controls. May be empty.}
  \item{\code{bin.pars}:}{Integer \code{vector} indicating
  which parameters to use for binning. These must be in the
  same position in all tubes.}
  \item{\code{measure.pars}:}{\code{list} of integer
  \code{vector}s indicating which parameters to use for
  measurement. These must be specified per tube.} Note: all
  slots can be get and set using accessor methods, for
  example \code{bin.pars(myFlowSet) <- c(1,2,5)} }
}

