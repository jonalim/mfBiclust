# Use holdoutRepeat = NULL to do as many folds as max(dim(m) ./ holdouts) aka
# k-fold but with random (not indexed) fold sampling each time. Known from esaBcv v1.2.1 source code.
rcvs <- function(m, maxPCs, holdoutRepeat = 100, niter = 1, center = FALSE) {
  # EsaBcv(r.limit = 20, nRepeat = 12, center = FALSE)
  eb <- esaBcv::EsaBcv(Y = m, niter = niter, r.limit = maxPCs, 
                       nRepeat = holdoutRepeat, center = center, only.r = TRUE)
  eb$result.list
}

bcv <- function(m, maxPCs, repetition = 20, center = TRUE) {
  # ps <- unlist(apply(rcvs, MARGIN = 1, function(row) {
  #   
  #   min(which(row < 0)) - 1
  #   which.min(row)
  # }))
  
  ebTest <- esaBcv::EsaBcv(Y = m, r.limit = 3, nRepeat = NULL)
  k <- min(nrow(ebTest$result.list), ncol(m), 10)
  
  ps <- unlist(sapply(seq_len(repetition), function(i) {
    # takes the argmin[index] of the colmeans of BCV prediction error
    results <- rcvs(m, maxPCs, holdoutRepeat = k^2, center = center)
    as.numeric(names(which.min(colMeans(results))))
  }))
  ps
}

evalEsaBcv.sim <- function(numBiclust = NULL, maxPCs = 10, center = TRUE, 
                           noise = 0.25, save = FALSE) {
  if (!is.null(numBiclust)) {
    sim <- genSimData(numBiclust, noise)
  } else {
    numBiclust <- "3_orig"
    sim <- genSimData3()
  }

  evalEsaBcv.matrix(sim, center, maxPCs, numBiclust, save)
}

evalEsaBcv.matrix <- function(m, center = FALSE, maxPCs = 10, fileid = "", 
                              save = FALSE) {

  if (save) { png(paste0("clusters", fileid, ".png"))}
  #plot
  old.par <- par(no.readonly = T)
  par(mar = c(0, 0, 0, 0))
  image(m, useRaster = TRUE, axes = FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), 
         c(min(m), round(max(m), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  if (save) { dev.off() }

  niter <-  1
  ebTest <- esaBcv::EsaBcv(Y = m, center = center, niter = niter, r.limit = 2, 
                           nRepeat = NULL)
  k <- min(nrow(ebTest$result.list), ncol(m), 10)
  res <- rcvs(m = m, maxPCs = maxPCs, holdoutRepeat = k ^ 2, center = center, 
              niter = niter)

  if (save) { png(paste0("bcvPE", fileid, ".png")) }
  df <- data.frame(bcv.PredictionError = rep(as.numeric(colnames(res)), 
                                             each = nrow(res)),
                   k = as.vector(res))
  boxplot(k ~ bcv.PredictionError, data = df)
  if (save) { dev.off() }
  
  as.numeric(names(which.min(colMeans(res))))
}

evalEsaBcv.file <- function(center = TRUE, save = TRUE) {
  data <- chooseFile()
  evalEsaBcv.matrix(data, center, 20, "_simdata5", save)
}


testCenter <- function(input) {
  genData3 <- genSimData(3, 0.01)
  simData3 <- genSimData3()

  uncentered <- sapply(1:1000, FUN = function(x) {
    sd3 <- tryCatch(
      esaBcv::EsaBcv(simData3, center = FALSE, niter = 1, r.limit = 6, 
                     nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    gd3 <- tryCatch(
      esaBcv::EsaBcv(genData3, center = FALSE, niter = 1, r.limit = 6, 
                     nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    i <- tryCatch(
      esaBcv::EsaBcv(input, center = FALSE, niter = 1, r.limit = 6, 
                     nRepeat = 2, only.r = TRUE)$best.r,
                  error = function(e) "error")
    c(sd3, gd3, i)
  })

  centered <- sapply(1:1000, FUN = function(x) {
    sd3 <- tryCatch(
      esaBcv::EsaBcv(simData3, center = TRUE, niter = 1, r.limit = 6, 
                     nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    gd3 <- tryCatch(
      esaBcv::EsaBcv(simData3, center = TRUE, niter = 1, r.limit = 6, 
                     nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    i <- tryCatch(
      esaBcv::EsaBcv(simData3, center = TRUE, niter = 1, r.limit = 6, 
                     nRepeat = 2, only.r = TRUE)$best.r,
                  error = function(e) "error")
    c(sd3, gd3, i)
  })

  uncentered <- t(uncentered)
  colnames(uncentered) <- c("simdata3", "gendata3", "simdata5")
  apply(uncentered, MARGIN = 2, table)

  centered <- t(centered)
  colnames(centered) <- c("simdata3", "gendata3", "simdata5")
  apply(centered, MARGIN = 2, table)
}