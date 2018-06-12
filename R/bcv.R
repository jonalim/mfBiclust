# Use holdoutRepeat = NULL to do as many folds as max(dim(m) ./ holdouts) aka
# k-fold but with random (not indexed) fold sampling each time. Known from esaBcv v1.2.1 source code.
rcvs <- function(m, maxPCs = 5, holdoutRepeat = 100, niter = 1, center = FALSE) {
   
}

esabcvWrapper <- function(m, maxPCs, holdoutRep = 20, niter = 1, center = TRUE) {
  res <- unlist(sapply(seq_len(repetition), function(x) {
    # takes the argmin[index] of the colmeans of BCV prediction error
    eb <- NULL
    while(is.null(eb)) {
      try(
        eb <- esaBcv::EsaBcv(Y = m, niter = niter, r.limit = maxPCs, nRepeat = holdoutRep, center = center, only.r = TRUE)
      )
    }
    
    as.numeric(names(which.min(eb$result.list[1, ])))
  }))
  
  list(best = mean(res), results = res)
}

evalEsaBcv.sim <- function(numBiclust = NULL, maxPCs = 10, center = TRUE, noise = 0.25, save = FALSE) {
  if(!is.null(numBiclust)) { 
    sim <- genSimData(numBiclust, noise)
  } else {
    numBiclust <- "3_orig"
    sim <- genSimData3()
  }
  
  evalEsaBcv.matrix(sim, center, maxPCs, numBiclust, save)
}

evalEsaBcv.matrix <- function(m, center = FALSE, maxPCs = 10, fileid = "", save = FALSE) {
  
  if(save) { png(paste0("clusters", fileid, ".png"))}
  #plot
  old.par <- par(no.readonly=T) 
  par(mar=c(0, 0, 0, 0))
  image(m, useRaster=TRUE, axes=FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), c(min(m), round(max(m), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  if(save) {dev.off()}
  
  niter = 1
  ebTest <- esaBcv::EsaBcv(Y = m, center = center, niter = niter, r.limit = 2, nRepeat = NULL)
  k <- min(nrow(ebTest$result.list), ncol(m), 10)
  res <- rcvs(m = m, maxPCs = maxPCs, holdoutRepeat = k^2, center = center, niter = niter)
  
  if(save) { png(paste0("bcvPE", fileid, ".png")) }
  df <- data.frame(bcv.PredictionError = rep(as.numeric(colnames(res)), each = nrow(res)),
                   k = as.vector(res))
  boxplot(k ~ bcv.PredictionError, data = df)
  if(save) {dev.off()}
  
  as.numeric(names(which.min(colMeans(res))))
}

evalEsaBcv.file <- function(center = TRUE, save = TRUE) {
  data <- chooseFile()
  # when I performed this 20 times on simdata5, the results were
  # 16, 14, 11, 12, 9, 11, 13 and then crash
  # I think it is a bug in esaBcv that when the variance condition isn't satisfied
  # on the k = 2, the whole function crashes
  evalEsaBcv.matrix(data, center, 20, "_simdata5", save)
}


testCenter <- function(input) {
  genData3 <- genSimData(3,0.01)
  simData3 <- genSimData3()
  
  iter1 <- sapply(1:1000, FUN = function(x) {
    sd3 <- tryCatch(esaBcv::EsaBcv(simData3, center = FALSE, niter = 1, r.limit = 6, nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    gd3 <- tryCatch(esaBcv::EsaBcv(genData3, center = FALSE, niter = 1, r.limit = 6, nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    c(sd3, gd3)
  })
  
  iter3 <- sapply(1:1000, FUN = function(x) {
    sd3 <- tryCatch(esaBcv::EsaBcv(simData3, center = FALSE, niter = 3, r.limit = 6, nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    gd3 <- tryCatch(esaBcv::EsaBcv(simData3, center = FALSE, niter = 3, r.limit = 6, nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    c(sd3, gd3)
  })
  
  iter1 <- t(iter1)
  colnames(iter1) <- c("simdata3", "gendata3")
  apply(iter1, MARGIN = 2, table)
  
  iter3 <- t(iter3)
  colnames(iter3) <- c("simdata3", "gendata3")
  apply(iter3, MARGIN = 2, table)
}