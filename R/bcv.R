# Use holdoutRepeat = NULL to do as many folds as max(dim(m) ./ holdouts) aka
# k-fold but with random (not indexed) fold sampling each time. Known from esaBcv v1.2.1 source code.
rcvs <- function(m, maxPCs = 5, holdoutRepeat = 100, niter = 1, center = FALSE) {
  
}

esabcvWrapper <- function(m, maxPCs, holdoutRep = 2, niter = 3, center = TRUE) {
  res <- unlist(sapply(seq_len(holdoutRep), function(x) {
    # takes the argmin[index] of the colmeans of BCV prediction error
    eb <- NULL
    
    while(is.null(eb)) {
      try(
        eb <- esaBcv::EsaBcv(Y = m, niter = niter, r.limit = maxPCs, nRepeat = 2, center = center, only.r = TRUE)
      )
    }
    
    as.numeric(names(which.min(eb$result.list[1, ])))
  }))
  
  list(best = mean(res), results = res)
}

analyzeTest <- function(test_res, k_limit = 10, iter = 10) {
  lapply(test_res, FUN = function(l) {
    esaBcvRes <- l[[1]]
    as.numeric(names(which.min(colMeans(esaBcvRes, na.rm = TRUE))))
    try(plot(rep(0:k_limit, times = iter), as.vector(esaBcvRes)))
    
    browser()
    
    bcvRes <- l[[2]]
    mean(bcvRes)
    plot(rep(1, iter), bcvRes)
    browser()

  })
  
}

testBcvNewOld <- function(directory, k_limit = 10, iter = 10) {
  res <- lapply(list.files(directory), function(file) {
    path <- file.path(directory, file)
    
    input <- read.csv(file = path)
    esaRes <- sapply(1:iter, FUN = function(x) {
      
      esabcv <- tryCatch({
        eb <- esaBcv::EsaBcv(input, center = TRUE, niter = 3, r.limit = k_limit, nRepeat = 2, only.r = TRUE)$result.list[1, ]
        if(length(eb) < k_limit + 1) { 
          eb <- c(eb, rep(0, each = k_limit + 1 - length(eb))) 
        }
        eb
      },
      error = function(e) rep(NA, each = 11))
      esabcv
    })
    esaRes <- t(esaRes)
    # returns table of test statistic itself
    
    bcvRes <- sapply(1:iter, FUN = function(x) {
      which.min(bcv(input, k_limit))
    })
    # returns just a vector of the results
    res <- list(esaRes, bcvRes)
    names(res) <- c("eb", "b")
    res
  })
  names(res) <- list.files(directory)
  res
}
# 
# eval_bcv <- function(Y, k_limit, repeats = 20, holdouts = 10) {
#   results <- sapply(seq_len(repeats), function(x) {bcv(simdata3, k_limit = k_limit)})
#   k.best <- mean(apply(results, MARGIN = 2, FUN = function(col) { which.min(col) }))
#   list(k.best = k.best, results = results)
# }

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
