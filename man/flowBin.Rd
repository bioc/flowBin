\name{flowBin}
\alias{flowBin}
\alias{flowbin-package}
\title{function to run the entire flowBin pipeline}
\usage{
  flowBin(tube.list, bin.pars, control.tubes = vector(),
    measure.pars = NULL,
    sample.name = "Unnamed Flow Expr Set",
    bin.method = "kmeans", expr.method = "medianFI",
    sparse.bin.thresh = 0.001, dequantize = T,
    snow.cluster = NULL, n.bins = 128, scale.expr = F,
    do.qnorm = T, return.bins = F)
}
\arguments{
  \item{tube.list}{a list of flowFrames, one for each tube
  to combine}

  \item{bin.pars}{a numerical vector indicating which flow
  parameters in the each flowFrame to use for combining
  tubes. These should be the same markers assayed across
  every tube.}

  \item{control.tubes}{a vector indicating which tubes in
  tube.list to use for negative controls. May be empty.}

  \item{measure.pars}{a list of which parameters to measure
  expression for, with one vector for each tube. If left
  NULL, this defaults to all parameters other than those
  specified as \code{bin.pars}}

  \item{sample.name}{name of this flowSample, for
  convenience (defaults to 'Unnamed Flow Expr Set')}

  \item{bin.method}{The method to use for creating bins.
  The two options are "kmeans" for k-means clustering and
  nearest-neighbour mapping of bins. or "flowFP" for flowFP
  binning and direct mapping of bin boundaries across
  tubes.}

  \item{expr.method}{The method to use to compute bin
  expression across tubes. This defaults to MFI of the
  cells belonging to that bin in each tube. Other options
  are}

  \item{sparse.bin.thresh}{Bins which contain fewer than
  this proportion of total events in any tube will be
  excluded as outliers. Defaults to 0.001}

  \item{dequantize}{If TRUE, adds a small (region of 1e-8)
  value to flow data to help break ties when binning.}

  \item{snow.cluster}{A cluster created using the
  \code{snow} package, which flowBin will use to speed up
  computation. If NULL, flowBin will execute in serial
  mode.}

  \item{n.bins}{Number of bins to use. Note that this must
  be a power of 2 if flowFP is selected as binning method.}

  \item{scale.expr}{If TRUE, the resulting expression
  values will be scaled to (0,1) using the ranges specified
  in the flowFrames in \code{tube.list}.}

  \item{do.qnorm}{If TRUE, the binning markers will be
  quantile normalized prior to binning.}

  \item{return.bins}{If TRUE, return a
  \code{BinnedFlowExprSet} containing the bins themselves
  as well as the expression for each bin.}
}
\value{
  A matrix containing expression values for each bin in
  terms of each marker across all tubes. If
  \code{return.bins} is set TRUE, then a list containing a
  \code{BinnedFlowExprSet} followed by the expression
  matrix is returned.
}
\description{
  Takes a list of flowFrames representing tubes from a
  single flow cytometry sample, and combines them using
  binning of events in terms of common markers across
  tubes.
}
\examples{
data(amlsample)
tube.combined <- flowBin(aml.sample@tube.set,
bin.pars=aml.sample@bin.pars, 
bin.method='flowFP',
control.tubes=aml.sample@control.tubes, 
expr.method='medianFIDist', scale.expr=TRUE)
heatmap(tube.combined, scale='none')
}

