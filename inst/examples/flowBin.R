data(amlsample)
tube.combined <- flowBin(aml.sample@tube.set,
bin.pars=aml.sample@bin.pars, 
bin.method='flowFP',
control.tubes=aml.sample@control.tubes, 
expr.method='medianFIDist', scale.expr=TRUE)
heatmap(tube.combined, scale='none')
