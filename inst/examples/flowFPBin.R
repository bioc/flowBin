data(amlsample)
normed.sample <- quantileNormalise(aml.sample)
res <- flowFPBin(normed.sample)