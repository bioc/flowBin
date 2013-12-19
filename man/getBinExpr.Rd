\name{getBinExpr}
\alias{getBinExpr}
\alias{getBinExpr,BinnedFlowSample-method}
\title{Calculate the expression of each bin in a BinnedFlowSample in terms of the measurement markers}
\arguments{
  \item{method}{Method to use to compute expression, passed
  as a string. Defaults to \code{medianFI}, which takes the
  simple median of each bin, and does not require control
  tubes. Other options available are \code{medianFIDist},
  which uses medians with the median of the negative
  control subtracted out, and \code{propPos} which sets a
  threshold at the 98th percentile of the negative control
  and determines what proportion of cells lie above that.}

  \item{include.bin.medians}{logical, specifies whether to
  compute the medians of each bin in terms of the binning
  markers and include them in the result or not. Defaults
  to \code{T}.}

  \item{scale}{logical specifying whether to scale the
  results to the interval (0,1). If \code{T} (default),
  then all medians will be divided by the range for that
  marker as specified in the flowFrame.}
}
\value{
  A numeric matrix containing expression values, with bins
  as rows and markers as columns
}
\description{
  getBinExpr main function definition
}
\examples{
data(amlsample)
normed.sample <- quantileNormalise(aml.sample)
binned.sample <- flowFPBin(normed.sample)
binned.sample <- removeSparseBins(binned.sample, 0.001)
bin.expr <- getBinExpr(binned.sample)
heatmap(bin.expr, scale='none')
}

