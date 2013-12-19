\name{removeSparseBins}
\alias{removeSparseBins}
\alias{removeSparseBins,BinnedFlowSample-method}
\title{Remove bins from a BinnedFlowSample with few events in them}
\arguments{
  \item{object}{the BinnedFlowSample to act on}

  \item{cutoff.prop=NULL}{the minimum proportion that a bin
  must contain to be kept. If NULL, only bins with no
  events in at least one tube will be removed.}
}
\value{
  a BinnedFlowSample with sparse bins removed
}
\description{
  Remove bins from a BinnedFlowSample with few events in
  them
}
\details{
  This is important to do prior to calculating bin
  expression, as bins containing 2 or less events, for,
  example, cannot have their median computed.
}

