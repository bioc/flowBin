\name{flowFPBin}
\alias{flowFPBin}
\alias{flowFPBin,FlowSample-method}
\title{Bin sample using flowFP binning}
\arguments{
  \item{object}{flowSample to bin}

  \item{n.bins=128}{number of bins to use. This should be a
  power of 2, and will be rounded down to the nearest power
  of 2 if not.}

  \item{snow.cluster=NULL}{Optional snow cluster to use for
  parallel execution.}

  \item{dequantize=T}{If TRUE, adds a small (region of
  1e-8) value to flow data to help break ties when
  binning.}
}
\value{
  a \code{BinnedFlowSample}
}
\description{
  Bin sample using flowFP binning
}
\examples{
data(amlsample)
normed.sample <- quantileNormalise(aml.sample)
res <- flowFPBin(normed.sample)
}

