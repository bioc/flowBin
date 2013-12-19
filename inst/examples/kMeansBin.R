data(amlsample)
normed.sample <- quantileNormalise(aml.sample)
res <- kMeansBin(normed.sample)