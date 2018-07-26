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
#' @param ks a vector of bicluster quantities to consider
#' @param maxIter maximum number of iterations
#' @param tol tolerance used to determine convergence
#' @param bestOnly if FALSE, both the predicted number of biclusters and a table
#'   of result counts is returned
#'
#' @export
auto_bcv <- function(Y, ks, holdouts = 10L, maxIter = 100L, tol = (10 ^ -4), bestOnly = TRUE,
                     verbose = TRUE, duplicable = TRUE, interactive = TRUE) {
  oldSeed <- duplicable("autobc") # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  # set up variables for testing convergence of the results
  distr <- rep(1L, each = length(ks))
  names(distr) <- as.character(ks)
  distrOld <- distr
  i <- 0L
  converged <- FALSE
  while(!converged && i < maxIter) {
    # Get the number of biclusters with lowest bcv value
    res <- bcv(Y, ks, duplicable = FALSE, holdouts = holdouts, interactive)
    bcvRes <- names(which.min(res))
    distr[bcvRes] <- distr[bcvRes] + 1L
    
    # In case some of the highest ks were not tested by bcv, amend distr.
    # This should occur on the first iteration only.
    if(length(distr) > length(res)) distr <- distr[seq_along(res)]
    names(distr) <- names(res)
    
    resid <- unlist(mapply(function(d, dOld) {
      (d / sum(distr) - dOld / sum(distrOld)) ^ 2
    }, d = distr, dOld = distrOld))
    if(verbose) {
      cat(paste("Iteration", i + 1L, "\n")) 
      cat("BCV result distribution:\n")
      message(paste0(
        do.call(paste, as.list(c(as.list(names(distr)), sep = "\t"))), "\n",
        do.call(paste, as.list(c(as.list(distr - 1L), sep = "\t")))
      ))
    }
    converged <- all(resid < tol)
    
    distrOld <- distr
    i <- i + 1L
  }
  if(i == maxIter) {
    warning("BCV results did not converge after", maxIter, "iterations")
  }
  
  distr <- distr - 1 # remove the pseudocount that was added
  
  med <- names(distr[min(which(cumsum(distr) > (sum(distr) / 2)))])
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
#' biclusters evaluated. The range of biclusters may be non-contiguous.
#'
#' @param Y the input matrix
#' @param ks the range of biclusters to evaluate
#' @param holdouts the number of holdouts to perform along each dimension of
#'   matrix \code{Y}
#'
#' @export
bcv <- function(Y, ks, holdouts = 10L, duplicable = TRUE, interactive = TRUE) {
  if(duplicable) {
    oldSeed <- duplicable("bcv") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  }
  
  #validate input
  if(!inherits(Y, "matrix")) Y <- as.matrix(Y)
  if(any(is.na(Y)) && interactive) {
    message(paste("Since some elements of the matrix are NA, BCV will run the",
                  "slower iterative algorithm. Continue?"))
    if (menu(c("Yes", "No")) != 1) {
      stop("User choice")
    }
  }
  p <- ncol(Y)
  n <- nrow(Y)
  
  kLimit <- floor((holdouts - 1L) / holdouts * min(p, n))
  if(sum(ks < kLimit) < length(ks)) {
    warning(paste("Due to holdouts, some of the highest requested ks will not be",
                  "tested."))
  }
  ks <- ks[ks < kLimit]
  if(length(ks) < 1) { stop(paste("ks must be a range of integers less than the",
                                  "smaller matrix dimension"))
  }
  
  # Automatically decrease the number of holdouts if necessary
  # Ex: for dimension = 8, holdouts = 10, the result is 8
  if(holdouts > min(p, n)) {
    holdouts <- min(p, n)
    warning(paste("Using", holdouts, "holdouts to ensure non-zero holdout",
                  "size."))
  }
  
  result.list <- bcvGivenKs(Y, ks, holdouts)
  names(result.list) <- as.character(ks)
  
  result.list
}

###% Algorithm by A.B. Owen, P.O. Perry. Bi-cross-validation of the SVD and the 
###% nonnegative matrix factorization. Ann. Appl. Stat., 3 (2009), pp. 564-594
# Helper function
bcvGivenKs <- function(Y, ks, holdouts = 10L) {
  # initialize...
  p <- ncol(Y)
  n <- nrow(Y)
  
  # Set up bi-cross-validation folds
  nHoldoutInd <- (seq_len(n) %% holdouts) + 1L
  nHoldoutInd <- sample(nHoldoutInd, size = n)
  
  pHoldoutInd <- (seq_len(p) %% holdouts) + 1L
  pHoldoutInd <- sample(pHoldoutInd, size = p)
  
  # Try NIPALS-PCA if any NA values
  if(any(is.na(Y))) {
    warning(paste("Using NIPALS-PCA because some matrix elements are NA",
                  "This feature might fail if too many elements are NA."))
    pca <- function(Y, k) { 
      res <- nipals_pca_autoclean(Y, k, center = TRUE)$genericFit
      list(scores = res@fit@W, loadings = res@fit@H)
    }
  } else {
    pca <- function(Y, k) {
      res <- prcomp(Y, rank. = k, retx = TRUE, center = TRUE)
      list(scores = res$x, loadings = t(res$rotation))
    }
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
        estD_k <- MASS::ginv(tcv[, 1L:k, drop = FALSE] %*% pcv[1L:k, , drop = FALSE])
        
        # Approximation of the hold-out quadrant
        estA <- Y[rInd, -sInd, drop = FALSE] %*% estD_k %*% Y[-rInd, sInd, drop = FALSE]
        
        resid <- A - estA # residual holdout matrix
        sum(resid ^ 2, na.rm = TRUE) # squared Frobenius norm
      })
      holdoutRes
    })
    
    rowSums(colHoldoutResults)
  })
  rcvs <- rowSums(rowHoldoutResults)
  names(rcvs) <- as.character(ks)
  
  # In all, we took the sum of residuals of A - estA, for holdouts covering A
  # exactly once. Therefore, we normalize to the size of the matrix here,
  # without NA's because NAs are evaluated as 0 when taking the Frobenius norm.
  rcvVals <- rcvs / (p * n - sum(is.na(Y)))
  
  # subtract <- c(0, cumsum(var_k[1:(length(var_k) - 1)]))
  # frob_square_x <- sum(Y ^ 2)
  # rcvVals <- mapply(function(rcv, subtract, var_k) {
  #   rcv <- 100 * (1 - rcv / frob_square_x)
  #   rcv <- rcv - subtract
  #   min(rcv, var_k)
  # }, rcv = rcvs, subtract = subtract, var_k = var_k, SIMPLIFY = TRUE)
  rcvVals
}
