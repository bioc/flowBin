\name{kMeansBin}
\alias{kMeansBin}
\alias{kMeansBin,FlowSample-method}
\title{Bin sample using K-means binning}
\arguments{
  \item{object}{flowSample to bin}

  \item{n.bins=128}{number of bins to use. This should be a
  power of 2, and will be rounded down to the nearest power
  of 2 if not.}

  \item{n.neighbours=1}{number of neighbours to use for KNN
  mapping of bins from clustered tube}

  \item{snow.cluster=NULL}{Optional snow cluster to use for
  parallel execution.}

  \item{random.seed=101}{Random seed to set to make K-means
  clustering deterministic.}

  \item{dequantize=T}{If TRUE, adds a small (region of
  1e-8) value to flow data to help break ties when
  binning.}
}
\value{
  a \code{BinnedFlowSample}
}
\description{
  Bin sample using K-means binning
}
\details{
  Runs K-means clustering on the binning markers in the
  first tube of the data set. These clusters are then
  mapped to the other tubes using K-nearest neighbours.
}
\examples{
data(amlsample)
normed.sample <- quantileNormalise(aml.sample)
res <- kMeansBin(normed.sample)
}

