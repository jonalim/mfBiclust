#' Use BCV to estimate the optimal k until convergence
#'
#' Repeatedly chooses the k that minimizes the BCV until the the distribution of
#' k's has converged. Often less than 50 iterations are required. The median k
#' is reported as "best".
#'
#' The highest k tested is limited to \eqn{(\code{holdouts} - 1) /
#' \code{holdouts} * \min(m, n)} for \eqn{Y_{m,n}}. A warning will be issued if
#' not all \code{ks} can be tested.
#'
#' @param Y the input matrix
#' @param ks a vector of bicluster quantities to consider
#' @param holdouts the number of row and column partitions. The true number of
#'   holdouts will be \code{holdouts} ^ 2.
#' @param maxIter maximum number of iterations
#' @param tol tolerance used to determine convergence
#' @param bestOnly if FALSE, both the predicted number of biclusters and a table
#'   of result counts is returned
#' @param verbose provide output after each iteration
#' @param duplicable fix the random seed internally
#' @param interactive prompt before running bcv on matrices with missing values
#'
#' @return if \code{bestOnly = FALSE}, the predicted bicluster quantity. if
#'   \code{bestOnly = TRUE}, a \code{\link{list}} containing: \describe{
#'   \item{best}{the predicted bicluster quantity} \item{counts}{a named table
#'   of result counts} }
#' @seealso \code{\link{bcv}()}
#' 
#' @examples
#' auto_bcv(yeast_benchmark[[1]])
#' @export
auto_bcv <- function(Y, ks = seq_len(min(nrow(Y), ncol(Y)) - 1), holdouts = 3L,
                     maxIter = 100L, tol = (10 ^ -4), bestOnly = FALSE,
                     verbose = FALSE, duplicable = TRUE, interactive = TRUE) {
    
    # Would be nice to make this a generic that can be called on matrix or
    # BiclusterExperiment objects

    oldSeed <- duplicable("autobc") # do not modify the R global environment
    on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    
    #validate input
    if(!inherits(Y, "matrix")) Y <- as.matrix(Y)
    if(any(is.na(Y)) && interactive) {
        message(paste("Since some elements of the matrix are NA, BCV will run",
                      "a much slower iterative algorithm. Continue?"))
        if (menu(c("Yes", "No")) != 1) {
            stop("User choice")
        }
    }
    # validate ks with warning if any are invalid
    ks <- kLimiter(holdouts, min(ncol(Y), nrow(Y)), ks)
    
    # set up variables for testing convergence of the results
    distr <- rep(1L, each = length(ks))
    names(distr) <- as.character(ks)
    distrOld <- distr
    i <- 0L
    converged <- FALSE
    while(!converged && i < maxIter) {
        # Get the number of biclusters with lowest bcv value
        res <- bcv(Y, ks, duplicable = FALSE, holdouts = holdouts,
                   interactive = FALSE)
        bcvRes <- names(which.min(res))
        distr[bcvRes] <- distr[bcvRes] + 1L
        
        # In case some of the highest ks were not tested by bcv, amend distr.
        # This should occur on the first iteration only.
        if(length(distr) > length(res)) distr <- distr[seq_along(res)]
        names(distr) <- names(res)
        distrOld <- distrOld[seq_along(distr)] # in case any ks were rejected
        
        resid <- unlist(mapply(function(d, dOld) {
            (d / sum(distr) - dOld / sum(distrOld)) ^ 2
        }, d = distr, dOld = distrOld))
        if(verbose) {
            cat(paste("Iteration", i + 1L, "\n"))
            cat(paste("Current best:",
                      names(which.max(distr)),
                      "\n"))
            cat("BCV result distribution:\n")
            message(paste0(
                do.call(paste, as.list(c(as.list(names(distr)), sep = "\t"))),
                "\n",
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
    
    med <- as.numeric(names(which.max(distr)))
    if(bestOnly) { med }
    else { list(best = med, counts = distr) }
}

#' Calculate bi-cross-validation
#'
#' Calculates the bi-cross-validation of a Singular value decmoposition of
#' matrix \code{Y}.
#'
#' It is recommended to use auto_bcv to perform and analyze several replications
#' of bi-cross-validation.
#'
#' @param Y the input matrix
#' @param ks a vector of bicluster quantities to evaluate
#' @param holdouts the number of row and column partitions. The true number of
#'   holdouts will be \code{holdouts} ^ 2.
#' @param duplicable fix the random seed internally
#' @param interactive prompt before running if \code{Y} is missing values
#'
#' @return A named vector of BCV values corresponding to the various numbers of
#'   biclusters evaluated.
#' @seealso \code{\link{auto_bcv}()}
#' 
#' @examples
#' bcv(yeast_benchmark[[1]])
#' 
#' @export
bcv <- function(Y, ks = seq_len(min(nrow(Y), ncol(Y)) - 1), holdouts = 10L, 
                duplicable = FALSE, interactive = TRUE) {
    if(duplicable) {
        oldSeed <- duplicable("bcv") # do not modify the R global environment
        on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    }
    
    #validate input
    if(!inherits(Y, "matrix")) Y <- as.matrix(Y)
    if(any(is.na(Y)) && interactive) {
        message(paste("Since some elements of the matrix are NA, BCV will run",
                      "a much slower iterative algorithm. Continue?"))
        if (menu(c("Yes", "No")) != 1) {
            stop("User choice")
        }
    }
    minDim <- min(ncol(Y), nrow(Y))
    ks <- kLimiter(holdouts, minDim, ks)
    
    if(length(ks) > 1 || is.numeric(ks)) {
        if (any(ks < 1) || any(!is.wholenumber(ks))) {
            stop(paste("ks must be a vector of positive whole numbers"))
        }
    } else {
        stop(paste("ks must be a vector of positive whole numbers"))
    }
    
    # Automatically decrease the number of holdouts if necessary
    # Ex: for dimension = 8, holdouts = 10, the result is 8
    if(holdouts > minDim) {
        holdouts <- minDim
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
        pca <- function(Y, k) { 
            res <- nipals_pca(Y, k, center = TRUE, scale = FALSE)$genericFit
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
    rowHoldoutResults <- vapply(
        seq_len(holdouts), FUN.VALUE = numeric(length(ks)), FUN = function(x) {
            rInd <- which(nHoldoutInd == x) # row holdout indices
            
            # Returns an array, ks by holdouts. Need to take row sums.
            colHoldoutResults <- vapply(
                seq_len(holdouts), FUN.VALUE = numeric(length(ks)),
                FUN = function(x) {
                    sInd <- which(pHoldoutInd == x) # column holdout indices
                    
                    A <- Y[rInd, sInd] # "holdout quadrant"
                    D <- Y[-rInd, -sInd] # "holdin quadrant"
                    
                    # PCA to factor the holdin quadrant
                    pcares11 <- pca(D, max(ks))
                    tcv <- pcares11$scores
                    pcv <- pcares11$loadings
                    
                    # Returns k norms. Must sum these up to obtain rcvs
                    holdoutRes <- vapply(
                        ks, FUN.VALUE = numeric(1), FUN = function(k) {
                            # PCA-based approximation of the hold-in quadrant
                            estD_k <- tcv[, 1L:k, drop = FALSE] %*%
                                pcv[1L:k, , drop = FALSE]
                            
                            # Approximation of the hold-out quadrant
                            estA <- Y[rInd, -sInd, drop = FALSE] %*%
                                MASS::ginv(estD_k) %*%
                                Y[-rInd, sInd, drop = FALSE]
                            
                            resid <- A - estA # residual holdout matrix
                            # return squared Frobenius norm
                            return(sum(resid ^ 2, na.rm = TRUE))
                        })
                    return(holdoutRes)
                })
            return(rowSums(colHoldoutResults))
        })
    rcvs <- rowSums(rowHoldoutResults)
    names(rcvs) <- as.character(ks)
    
    # In all, we took the sum of residuals of A - estA, for holdouts covering A
    # exactly once.
    return(rcvs)
}

kLimiter <- function(holdouts, dim, ks) {
    kLimit <- floor((holdouts - 1L) / holdouts * dim)
    if(sum(ks < kLimit) < length(ks)) {
        warning(paste("Due to the number of holdouts, not all ks can be",
                      "tested. A higher number of holdouts can be specified to",
                      "test higher values of k, but",
                      "the author does not recommend doing so.\n"),
                call. = FALSE)
    }
    return(ks[ks < kLimit])
}
