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
# 
# evalEsaBcv.sim <- function(numBiclust = NULL, maxPCs = 10, center = TRUE, 
#                            noise = 0.25, save = FALSE) {
#   if (!is.null(numBiclust)) {
#     sim <- genSimData(numBiclust, noise)
#   } else {
#     numBiclust <- "3_orig"
#     sim <- genSimData3()
#   }
# 
#   evalEsaBcv.matrix(sim, center, maxPCs, numBiclust, save)
# }

# evalEsaBcv.matrix <- function(m, center = FALSE, maxPCs = 10, fileid = "", 
#                               save = FALSE) {
# 
#   if (save) { png(paste0("clusters", fileid, ".png"))}
#   #plot
#   old.par <- par(no.readonly = T)
#   par(mar = c(0, 0, 0, 0))
#   image(t(apply(m, 2, rev)), useRaster = TRUE, axes = FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
#   legend(grconvertX(0.5, "device"), grconvertY(1, "device"), 
#          c(min(m), round(max(m), digits = 3)),
#          fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
#   par(old.par)
#   if (save) { dev.off() }
# 
#   niter <-  1
#   ebTest <- esaBcv::EsaBcv(Y = m, center = center, niter = niter, r.limit = 2, 
#                            nRepeat = NULL)
#   k <- min(nrow(ebTest$result.list), ncol(m), 10)
#   res <- rcvs(m = m, maxPCs = maxPCs, holdoutRepeat = k ^ 2, center = center, 
#               niter = niter)
# 
#   if (save) { png(paste0("bcvPE", fileid, ".png")) }
#   df <- data.frame(bcv.PredictionError = rep(as.numeric(colnames(res)), 
#                                              each = nrow(res)),
#                    k = as.vector(res))
#   boxplot(k ~ bcv.PredictionError, data = df)
#   if (save) { dev.off() }
#   
#   as.numeric(names(which.min(colMeans(res))))
# }
# 
# evalEsaBcv.file <- function(center = TRUE, save = TRUE) {
#   data <- chooseFile()
#   evalEsaBcv.matrix(data, center, 20, "_simdata5", save)
# }


testBcvNewOld <- function(directory, k_limit = 10, iter = 10) {
  res <- lapply(list.files(directory), function(file) {
    path <- file.path(directory, file)
    
    input <- read.csv(file = path)
    res <- sapply(1:iter, FUN = function(x) {
      
      esabcv <- tryCatch(esaBcv::EsaBcv(input, center = TRUE, niter = 3, r.limit = k_limit, nRepeat = 2, only.r = TRUE)$best.r,
                         error = function(e) "error")
      mybcv <- tryCatch({
        res <- bcv(input, k_limit)
        which.min(res)
      },
      error = function(e) "error")
      c(esabcv, mybcv)
    })
    rownames(res) <- c("eb", "b")
    res
  })
  
  names(res) <- list.files(directory)
  res
}

eval_bcv <- function(Y, k_limit, repeats = 20, holdouts = 10) {
  results <- sapply(seq_len(repeats), function(x) {bcv(simdata3, k_limit = k_limit)})
  k.best <- mean(apply(results, MARGIN = 2, FUN = function(col) { which.min(col) }))
  list(k.best = k.best, results = results)
}

bcv <- function(Y, k_limit, holdouts = 10) {
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
  
  nHoldoutInd <- sample(nHoldoutInd, size = n)
  pHoldoutInd <- sample(pHoldoutInd, size = p)
  
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
        # estD_kold <- tM[rInd, 1:k, drop = FALSE] %*% pM[1:k, sInd, drop = FALSE]
        estA <- Y[rInd, -sInd, drop = FALSE] %*% estD_k %*% Y[-rInd, sInd, drop = FALSE]
        
        # resid <- A - estD_kold - estA
        resid <- A - estA
        rcvs[x] <<- rcvs[x] + sum(resid ^ 2)
      })
    })
  })
  rcvVals <- rcvs / holdouts / holdouts # normalize for the number of holdouts
  
  # subtract <- c(0, cumsum(var_k[1:(length(var_k) - 1)]))
  # frob_square_x <- sum(Y ^ 2)
  # rcvVals <- mapply(function(rcv, subtract, var_k) {
  #   rcv <- 100 * (1 - rcv / frob_square_x)
  #   rcv <- rcv - subtract
  #   min(rcv, var_k)
  # }, rcv = rcvs, subtract = subtract, var_k = var_k, SIMPLIFY = TRUE)
  rcvVals
}
