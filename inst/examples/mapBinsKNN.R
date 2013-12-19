data(amlsample)
tube1.expr <- exprs(tube.set(aml.sample)[[1]])
kmeans.res <- kmeans(tube1.expr, 100)
kmeans.labels <- kmeans.res$cluster

#Now create a binnedFlowExprSet using the cluster labels for tube 1
clustered.sample <- mapBinsKNN(aml.sample, kmeans.labels)
sort(table(bin.labels(clustered.sample)[[3]]))