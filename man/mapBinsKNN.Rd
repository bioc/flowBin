\name{mapBinsKNN}
\alias{mapBinsKNN}
\alias{mapBinsKNN,FlowSample-method}
\title{Internal function to map bins by KNN}
\arguments{
  \item{object}{flowSample to map the bins of}

  \item{tube.1.labels}{integer vector of bin labels for the
  events in tube 1}

  \item{n.neighbours=1}{number of neighbours to use for KNN
  mapping of bins from clustered tube}

  \item{snow.cluster=NULL}{Optional snow cluster to use for
  parallel execution.}

  \item{dequant=T}{If TRUE, adds a small (region of 1e-8)
  value to flow data to help break ties when binning.}
}
\value{
  a \code{BinnedFlowSample}
}
\description{
  Internal function to map bins by KNN
}
\details{
  Takes a FlowSample and labels for the events in tube 1,
  and maps these to all other tubes.
}
\examples{
data(amlsample)
tube1.expr <- exprs(tube.set(aml.sample)[[1]])
kmeans.res <- kmeans(tube1.expr, 100)
kmeans.labels <- kmeans.res$cluster

#Now create a binnedFlowExprSet using the cluster labels for tube 1
clustered.sample <- mapBinsKNN(aml.sample, kmeans.labels)
sort(table(bin.labels(clustered.sample)[[3]]))
}

