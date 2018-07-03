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
auto_bcv <- function(Y, ks, maxIter = 100, tol = (10 ^ -4), bestOnly = TRUE,
                     verbose = TRUE) {
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  distr <- rep(1, each = length(ks))
  names(distr) <- as.character(ks)
  distrOld <- distr
  change <- 1
  i <- 0
  dold <- 0
  converged <- FALSE
  while(!converged && i < maxIter) {
    bcvRes <- which.min(bcv(Y, ks))
    distr[bcvRes] <- distr[bcvRes] + 1
    # all -> 340
    # sqrt(.Machine$double.eps) -> 230

    resid <- unlist(mapply(function(d, dOld) {
      (d / sum(distr) - dOld / sum(distrOld)) ^ 2
    }, d = distr, dOld = distrOld))
    if(verbose) print(resid)
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
bcv <- function(Y, ks, holdouts = 10) {
  oldSeed <- duplicable() # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  Y <- as.matrix(Y)
  p <- ncol(Y)
  n <- nrow(Y)
  
  ks <- ks[ks < min(nrow(Y) - 1, ncol(Y) - 1)]
  
  result.list <- bcvGivenKs(Y, ks, holdouts)
  names(result.list) <- as.character(ks)
  
  result.list
}

###% Algorithm by A.B. Owen, P.O. Perry. Bi-cross-validation of the SVD and the 
###% nonnegative matrix factorization. Ann. Appl. Stat., 3 (2009), pp. 564-594
# Helper function
bcvGivenKs <- function(Y, ks, holdouts = 10) {
  # initialize...
  p <- ncol(Y)
  n <- nrow(Y)

  # Set up bi-cross-validation folds
  nHoldoutInd <- (seq_len(n) %% holdouts) + 1
  nHoldoutInd <- sample(nHoldoutInd, size = n)
  
  pHoldoutInd <- (seq_len(p) %% holdouts) + 1
  pHoldoutInd <- sample(pHoldoutInd, size = p)

  # Try NIPALS-PCA if any NA values
  if(any(is.na(Y))) {
    pca <- function(Y, k) { nipals_pca(Y, k)$fit }
  } else {
    pca <- function(Y, k) { prcomp(Y, rank. = k, retx = TRUE, center = FALSE) }
  }
  
  # Perform PCA on the whole matrix
  pcares <- pca(Y, max(ks))
  tM = pcares$scores
  pM = pcares$loadings
  # var_k <- prcmp$sdev[seq_len(max(ks))] ^ 2
  
  # Returns an array, ks by holdouts. Need to take row sums again.
  rowHoldoutResults <- sapply(seq_len(holdouts), function(x) {
    rInd <- which(nHoldoutInd == x) # row holdout indices
    
    # Returns an array, ks by holdouts. Need to take row sums.
    colHoldoutResults <- sapply(seq_len(holdouts), function(x) {
      sInd <- which(pHoldoutInd == x) # column holdout indices
      
      A <- Y[rInd, sInd] # "holdout quadrant"
      D <- Y[-rInd, -sInd] # "holdin quadrant"
      
      # PCA to factor the holdin quadrant
      pcares11 <- pca(D, max(ks))
      tcv <- pcares11$scores
      pcv <- pcares11$loadings
      
      # Returns k norms. Must sum these up to obtain rcvs
      holdoutRes <- sapply(ks, function(k) {
        # PCA-based approximation of the hold-in quadrant
        estD_k <- MASS::ginv(tcv[, 1:k, drop = FALSE] %*% pcv[1:k, , drop = FALSE])
        
        # Approximation of the hold-out quadrant
        estA <- Y[rInd, -sInd, drop = FALSE] %*% estD_k %*% Y[-rInd, sInd, drop = FALSE]
        
        resid <- A - estA # residual holdout matrix
        sum(resid ^ 2) # squared Frobenius norm
      })
      holdoutRes
    })
    
    rowSums(colHoldoutResults)
  })
  rcvs <- rowSums(rowHoldoutResults)
  names(rcvs) <- as.character(ks)
  
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
