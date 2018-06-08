# Use holdoutRepeat = NULL to do as many folds as max(dim(m) ./ holdouts) aka
# k-fold but with random (not indexed) fold sampling each time. Known from esaBcv v1.2.1 source code.
rcvs <- function(m, maxPCs, holdoutRepeat = 100, center = FALSE) {
  # EsaBcv(r.limit = 20, nRepeat = 12, center = FALSE)
  eb <- esaBcv::EsaBcv(Y = m, r.limit = maxPCs, nRepeat = holdoutRepeat, center = center, only.r = TRUE)
  eb$result.list 
}

bcv <- function(m, maxPCs, repetition = 20, center = TRUE) {
  # ps <- unlist(apply(rcvs, MARGIN = 1, function(row) {
  #   
  #   min(which(row < 0)) - 1
  #   which.min(row)
  # }))
  
  ebTest <- esaBcv::EsaBcv(Y = m, r.limit = 3, nRepeat = NULL)
  k <- nrow(ebTest$result.list)
  
  ps <- unlist(sapply(seq_len(repetition), function(i) {
    # takes the argmin[index] of the colmeans of BCV prediction error
    results <- rcvs(m, maxPCs, holdoutRepeat = k^2, center = center)
    as.numeric(names(which.min(colMeans(results))))
  }))
  ps
}

evalEsaBcv.sim <- function(numBiclust, center = TRUE) {
  if(!is.null(numBiclust)) { 
    # genSimData for given number of biclusters
  } else {
    genData3 <- genSimData3()
  }
  
  Sys.sleep(5)
  
  res <- rcvs(m = genData3, maxPCs = 10, center = center)
  plot(x = rep(as.numeric(colnames(res)), each = nrow(res)), y = as.vector(res))

}

evalEsaBcv.file <- function(numBiclust = NULL, center = TRUE) {
  if(!is.null(numBiclust)) { 
    # genSimData for given number of biclusters
  }
  
  simdata5 <- chooseFile()
  
  res <- rcvs(m = genData3, maxPCs = 10, center = center)
  plot(x = rep(as.numeric(colnames(res)), each = nrow(res)), y = as.vector(res))
  
}