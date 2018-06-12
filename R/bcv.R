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

  iter1 <- sapply(1:1000, FUN = function(x) {
    sd3 <- tryCatch(esaBcv::EsaBcv(simData3, center = TRUE, niter = 1, r.limit = 6, nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    gd3 <- tryCatch(
      esaBcv::EsaBcv(genData3, center = TRUE, niter = 1, r.limit = 6, 
                     nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    c(sd3, gd3)
  })

  iter3 <- sapply(1:1000, FUN = function(x) {
    sd3 <- tryCatch(esaBcv::EsaBcv(simData3, center = TRUE, niter = 3, r.limit = 6, nRepeat = 2, only.r = TRUE)$best.r,
                    error = function(e) "error")
    gd3 <- tryCatch(esaBcv::EsaBcv(simData3, center = TRUE, niter = 3, r.limit = 6, nRepeat = 2, only.r = TRUE)$best.r,
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

bcv <- function(Y, k_limit, niter = 3, holdouts = 10) {
  Y <- as.matrix(Y)
  p <- ncol(Y)
  n <- nrow(Y)
  
  k_limit <- min(k_limit, nrow(Y) - 1, ncol(Y) - 1)
  
  result.list <- bcvGivenKs(Y, seq_len(k_limit), holdouts)
  names(result.list) <- seq_len(k_limit)
  
  result.list
}

bcvGivenKs <- function(Y, ks, holdouts = 10) {
  p <- ncol(Y)
  n <- nrow(Y)
  
  prcmp <- prcomp(Y, rank. = max(ks), retx = TRUE)
  tM = prcmp$x
  pM = t(prcmp$rotation)
  var_k <- prcmp$sdev[seq_len(max(ks))] ^ 2
  
  nHoldoutInd <- (seq_len(n) %% holdouts) + 1
  pHoldoutInd <- (seq_len(p) %% holdouts) + 1
  
  rcvs <- rep(0, length(ks))
  names(rcvs) <- as.character(ks)
  
  sapply(seq_len(holdouts), function(x) {
    # row holdout indices
    rInd <- which(nHoldoutInd == x) 
    
    sapply(seq_len(holdouts), function(x) {
      sInd <- which(pHoldoutInd == x)
      
      A <- Y[rInd, sInd]
      D <- Y[-rInd, -sInd] # 1/(holdouts^2) of the matrix X
      
      prcmp <- prcomp(D, rank. = max(ks), retx = TRUE)
      tcv <- prcmp$x
      pcv <- t(prcmp$rotation)
      
      
      sapply(seq_len(length(ks)), function(x) {
        k <- ks[x]
        estD_k <- MASS::ginv(tcv[, 1:k, drop = FALSE] %*% pcv[1:k, , drop = FALSE])
        estD_kold <- tM[rInd, 1:k, drop = FALSE] %*% pM[1:k, sInd, drop = FALSE]
        estA <- Y[rInd, -sInd, drop = FALSE] %*% estD_k %*% Y[-rInd, sInd, drop = FALSE]
        
        resid <- sum((A - estD_kold - estA) ^ 2)
        rcvs[x] <<- rcvs[x] + sum(resid ^ 2)
      })
    })
  })
  
  subtract <- c(0, cumsum(var_k[1:(length(var_k) - 1)]))
  frob_square_x <- sum(Y ^ 2)
  browser()
  rcvVals <- mapply(function(rcv, subtract, var_k) {
    rcv <- 100 * (1 - rcv / frob_square_x)
    rcv <- rcv - subtract
    min(rcv, var_k)
  }, rcv = rcvs, subtract = subtract, var_k = var_k, SIMPLIFY = TRUE)
  # FIXME rcvVals is really negative???
  rcvVals
}
