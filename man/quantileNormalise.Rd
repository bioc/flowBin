\name{quantileNormalise}
\alias{quantileNormalise}
\alias{quantileNormalise,FlowSample-method}
\title{quantileNormalise
normalise binning paramaters across all tubes of a flowSample}
\description{
  Since the binning parameters are the same across tubes,
  and samples each tube is an aliquot from the same sample,
  these should have the same underlying distribution.
  Hence, quantile normalisation can be used to force this
  to be so, removing technical variation.
}
\examples{
data(amlsample)
normed.sample <- quantileNormalise(aml.sample)
qnorm.check <- checkQNorm(aml.sample, normed.sample, do.plot=FALSE)
show(qnorm.check)
}

