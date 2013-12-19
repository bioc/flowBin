\name{checkQNorm}
\alias{checkQNorm}
\alias{checkQNorm,}
\alias{checkQNorm,FlowSample,FlowSample-method}
\title{Function to check the quantile normalisation of a FlowSample using flowFP}
\value{
  \code{list} containing two matrices of standard
  deviations across bins (rows) vs tubes (columns) for
  before (\code{sd.before}) and after (\code{sd.after}).
}
\description{
  Function to check the quantile normalisation of a
  FlowSample using flowFP
}
\details{
  object and normed.object are compared using \code{flowFP}
  binning, to assess the deviation in bin counts between
  the two.
}
\examples{
data(amlsample)
normed.sample <- quantileNormalise(aml.sample)
qnorm.check <- checkQNorm(aml.sample, normed.sample, do.plot=FALSE)
show(qnorm.check)
}

