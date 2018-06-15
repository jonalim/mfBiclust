# 
# 
# esabcvWrapper <- function(m, maxPCs, holdoutRep = 20, niter = 1, center = TRUE) {
#   res <- unlist(sapply(seq_len(repetition), function(x) {
#     # takes the argmin[index] of the colmeans of BCV prediction error
#     eb <- NULL
#     
#     while(is.null(eb)) {
#       try(
#         eb <- esaBcv::EsaBcv(Y = m, niter = niter, r.limit = maxPCs, nRepeat = holdoutRep, center = center, only.r = TRUE)
#       )
#     }
#     
#     as.numeric(names(which.min(eb$result.list[1, ])))
#   }))
#   
#   list(best = mean(res), results = res)
# }
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

# eval_bcv <- function(Y, k_limit, repeats = 20, holdouts = 10) {
#   results <- sapply(seq_len(repeats), function(x) {bcv(simdata3, k_limit = k_limit)})
#   k.best <- mean(apply(results, MARGIN = 2, FUN = function(col) { which.min(col) }))
#   list(k.best = k.best, results = results)
# }

#' Perform bcv until convergence
#'
#' Performs BCV until the the distribution of results has converged. Often this
#' requires less than 50 iterations.
#'
#' The returned number of biclusters is the median of the results from
#' \code{maxIter} iterations.
#'
#' @section Deciding the maximum k:
#' To tune kLimit, it might be helpful to run with \code{maxIter} around 10 and
#' \code{bestOnly = TRUE} to determine if there is an obvious upper bound on the
#' results.
#'
#' @param Y the input matrix
#' @param kLimit the maximum number of biclusters to consider
#' @param maxIter maximum number of iterations
#' @param tol tolerance used to determine convergence
#' @param bestOnly if FALSE, both the predicted number of biclusters and a table
#'   of result counts is returned
#'
#' @export
auto_bcv <- function(Y, kLimit, maxIter = 100, tol = (10 ^ -4), bestOnly = TRUE) {
  distr <- rep(1, each = kLimit)
  distrOld <- distr
  change <- 1
  i <- 0
  dold <- 0
  converged <- FALSE
  while(!converged && i < maxIter) {
    bcvRes <- which.min(bcv(Y, kLimit))
    distr[bcvRes] <- distr[bcvRes] + 1
    # all -> 340
    # sqrt(.Machine$double.eps) -> 230

    resid <- unlist(mapply(function(d, dOld) {
      (d / sum(distr) - dOld / sum(distrOld)) ^ 2
    }, d = distr, dOld = distrOld))
    #resid <- cumprod(resid)[length(resid)] ^ (1 / length(resid))
    print(resid)
    converged <- all(resid < tol)
    
    distrOld <- distr
    i <- i + 1
    
    if(i %% 10 == 0) { 
      print(paste("Iteration", i)) 
      print(paste("Highest residual was", max(resid)))
      }
  }
  if(i == maxIter) {
    warning("BCV results did not converge after", maxIter, "iterations")
  }
  
  med <- min(which(cumsum(distr) > (sum(distr) / 2)))
  if(bestOnly) { med }
  else { list(best = med, counts = distr) }
}
  
#' Perform bi-cross-validation
#'
#' The number of biclusters yielding the lowest BCV value is the predicted best.
#' It is recommended to use auto_bcv to perform several replications of
#' bi-cross-validation, and use the median of the results as the predicted
#' number of biclusters.
#'
#' A named vector of BCV values corresponding to the various numbers of
#' biclusters evaluated.
#'
#' @param Y the input matrix
#' @param kLimit the maximum number of biclusters to evaluate
#' @param holdouts the number of holdouts to perform along each dimension of
#'   matrix \code{Y}
#'
#' @export
bcv <- function(Y, kLimit, holdouts = 10) {
  Y <- as.matrix(Y)
  p <- ncol(Y)
  n <- nrow(Y)
  
  kLimit <- min(kLimit, nrow(Y) - 1, ncol(Y) - 1)
  
  result.list <- bcvGivenKs(Y, seq_len(kLimit), holdouts)
  names(result.list) <- seq_len(kLimit)
  
  result.list
}

# Helper function
bcvGivenKs <- function(Y, ks, holdouts = 10) {
  # initialize...
  p <- ncol(Y)
  n <- nrow(Y)
  
  rcvs <- rep(0, length(ks))
  names(rcvs) <- as.character(ks)
  
  # Set up bi-cross-validation folds
  nHoldoutInd <- (seq_len(n) %% holdouts) + 1
  nHoldoutInd <- sample(nHoldoutInd, size = n)
  
  pHoldoutInd <- (seq_len(p) %% holdouts) + 1
  pHoldoutInd <- sample(pHoldoutInd, size = p)

  # Perform PCA on the whole matrix
  prcmp <- prcomp(Y, rank. = max(ks), retx = TRUE)
  tM = prcmp$x
  pM = t(prcmp$rotation)
  var_k <- prcmp$sdev[seq_len(max(ks))] ^ 2
  
  sapply(seq_len(holdouts), function(x) {
    rInd <- which(nHoldoutInd == x) # row holdout indices
    
    sapply(seq_len(holdouts), function(x) {
      sInd <- which(pHoldoutInd == x) # column holdout indices
      
      A <- Y[rInd, sInd] # "holdout quadrant"
      D <- Y[-rInd, -sInd] # "holdin quadrant"
      
      # PCA to factor the holdin quadrant
      prcmp <- prcomp(D, rank. = max(ks), retx = TRUE)
      tcv <- prcmp$x
      pcv <- t(prcmp$rotation)
      
      sapply(seq_len(length(ks)), function(x) {
        k <- ks[x]
        # PCA-based approximation of the hold-in quadrant
        estD_k <- MASS::ginv(tcv[, 1:k, drop = FALSE] %*% pcv[1:k, , drop = FALSE])
        
        # Approximation of the hold-out quadrant
        estA <- Y[rInd, -sInd, drop = FALSE] %*% estD_k %*% Y[-rInd, sInd, drop = FALSE]
        
        resid <- A - estA # residual holdout matrix
        rcvs[x] <<- rcvs[x] + sum(resid ^ 2) # squared Frobenius norm
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
