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

evalEsaBcv.sim <- function(numBiclust = NULL, maxPCs = 10, center = TRUE, noise = 0.25, save = FALSE) {
  if(!is.null(numBiclust)) { 
    sim <- genSimData(numBiclust, noise, save = save)
  } else {
    numBiclust <- "3_orig"
    sim <- genSimData3(save = save)
  }
  
  evalEsaBcv.matrix(m, TRUE, maxPCs, numBiclust, save)
}

evalEsaBcv.matrix <- function(m, center = FALSE, maxPCs = maxPCs, fileid = "", save = FALSE) {
  
  if(save) { png(paste0("clusters", fileid, ".png"))}
  #plot
  old.par <- par(no.readonly=T) 
  par(mar=c(0, 0, 0, 0))
  image(m, useRaster=TRUE, axes=FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), c(min(m), round(max(m), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  if(save) {dev.off()}
  
  ebTest <- esaBcv::EsaBcv(Y = m, r.limit = 3, nRepeat = NULL, svd.method = "fast")
  k <- nrow(ebTest$result.list)
  
  res <- rcvs(m = m, maxPCs = maxPCs, holdoutRepeat = k^2, center = center)
  
  if(save) { png(paste0("bcvPE", fileid, ".png")) }
  df <- data.frame(bcv.PredictionError = rep(as.numeric(colnames(res)), each = nrow(res)),
                   k = as.vector(res))
  boxplot(k ~ bcv.PredictionError, data = df)
  if(save) {dev.off()}
  
  as.numeric(names(which.min(colMeans(res))))
}

evalEsaBcv.file <- function(center = TRUE, save = FALSE) {
  data <- chooseFile()
  
  evalEsaBcv.matrix(data, TRUE, 10, "_simdata5", save)
}