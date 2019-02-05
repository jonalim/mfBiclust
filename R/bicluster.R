#' Biclustering algorithms
#'
#' The \code{method} argument of \code{\link{addStrat}()} provides access to
#' internal \code{mfBiclust} functions, which themselves call functions from
#' packages \code{NMF} and \code{biclust} as back-ends.

#' Back-end arguments for "als-nmf", "svd-pca" and "nipals-pca" are described on
#' their respective pages. For arguments specific to other methods, please see
#' respective documentation in other packages.
#'
#' \describe{
#' \item{\link{als_nmf}}{Alternating-least-squares non-negative matrix
#' approximation. Fast at low values of \code{k}, but rapidly slows as \code{k}
#' increases.}
#' \item{\link{svd_pca}}{The Singular value decmoposition algorithm. Each
#' principal component is interpreted as the degree of membership in a single
#' bicluster. The resulting score matrix is thresholded to binarize bicluster
#' membership. \code{svd_pca} is the fastest provided algorithm.}
#' \item{\link{nipals_pca}}{An iterative PCA algorithm that may tolerate missing
#' data. Slower than \code{svd_pca}, but still faster than the other
#' algorithms.}
#' \item{plaid}{Uses \code{\link[=BCPlaid]{biclust::biclust(method = BCPlaid(),
#' ...)}} as its back-end. To get \code{k} biclusters, back-end arguments
#' \code{row.release} and \code{col.release} are simultaneously decremented
#' towards 0.1 in steps of 0.1, then \code{shuffle} is incremented towards 10 in
#' steps of 1.}
#' \item{snmf}{\strong{S}parse \strong{n}on-negative \strong{m}atrix
#' \strong{f}actorization. Uses
#' \code{\link[=nmfAlgorithm.SNMF_R]{NMF::nmf(method = "snmf/r", ...)}} as its
#' back-end. To get \code{k} biclusters, back-end arguments \code{eta} and
#' \code{beta} are initialized at mean(\code{A}) and are halved progressively.}
#' \item{spectral}{Uses \code{\link[=BCSpectral]{biclust::biclust(method =
#' BCSpectral(), ...)}} as its back-end. For smaller matrices, the back-end
#' argument \code{numberOfEigenvalues} is automatically set to enable finding
#' the number of biclusters requested. The back-end argument \code{withinVar} is
#' initialized equal to the smaller matrix dimension and allowed to increase up
#' to 10 times the smaller matrix dimension to find \code{k} biclusters.} }
#'
#' @seealso \code{\link{als_nmf}}
#' @seealso \code{\link{svd_pca}}
#' @seealso \code{\link[nipals]{nipals}}
#' @seealso \code{\link{nmfAlgorithm.SNMF_R}}
#' @seealso \code{\link[biclust]{BCSpectral}}
#' @seealso \code{\link[biclust]{BCPlaid}}
#' @name bicluster-methods
NULL

#' Alternating-least-squares non-negative matrix approximation
#'
#' Approximates a non-negative matrix as the product of two non-negative matrix
#' factors using the alternating-least-squares algorithm by Paatero and Tapper
#' (1994).
#'
#' Factorization is performed \code{reps} times, then the result with the
#' minimum mean squared-error is returned. \code{als_nmf()} is fast for a small
#' number of biclusters, but running time rapidly increases with \code{k}.
#'
#' @param A the matrix to factorize
#' @param k the number of factors to calculate
#' @param reps the number of replications to choose from
#' @param maxIter the maximum number of least-squares steps
#' @param eps_conv convergence tolerance
#' @param duplicable fix the random seed internally
#' @param verbose Print the mean squared error every 10 iterations
#' @param ... Present for compatibility. \code{als_nmf} has no other
#'   parameters
#'
#' @return a \code{\link{genericFit-class}} object
#' @export
#' @importFrom NMF .fcnnls
als_nmf <- function(A, k, reps = 4L, maxIter= 100L,
                    eps_conv = 1e-4, duplicable = TRUE, verbose= FALSE, ...){
    ###% Adapted from NMF v0.21.0 written by Renaud Gaujoux, Cathal Seoighe.
    ###% (2018)
    ###% https://cran.r-project.org/web/packages/NMF/
    ###% http://renozao.github.io/NMF
    if(duplicable) {
        oldSeed <- duplicable("biclus") # do not modify the R global environment
        on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    } else { reps <- 1 }

    if(any(A < 0)) {stop("Negative values are not allowed")}

    m = nrow(A); n = ncol(A); erravg1 = numeric();

    # to do with sparsity constraints
    maxA=max(A)

    ## VALIDITY of parameters
    # eps_conv
    if( eps_conv <= 0 )
        stop("ALS-NMF: Invalid argument 'eps_conv' - value should be positive")

    solutions <- lapply(seq_len(reps), function(i) {
        W <- NULL # these are the results of each replicate
        H <- NULL
        residNorm <- max(A)
        # initialize random W if no starting point is given

        # rank is given by k
        W <- matrix(runif(m*k), m,k)

        idxWold=rep(0, m); idxHold=rep(0, n); inc=0;

        # normalize columns of W
        frob <- apply(W, 2, function(x) sqrt(sum(x ^ 2)) )
        W <- sweep(W, 2, frob, FUN = "/")

        Wold <- W
        Hold <- matrix(runif(k*n), k,n);
        residNormOld <- maxA
        nrestart=0;
        restart <- function() {
            # re-initialize random W
            idxWold <<- rep(0, m); idxHold <<- rep(0, n); inc <<- 0;
            erravg1 <<- numeric();# re-initialize base average error
            W <- matrix(runif(m*k), m,k);

            # normalize columns of W
            frob <- apply(W, 2, function(x) sqrt(sum(x ^ 2)) )
            W <<- sweep(W, 2, frob, FUN = "/")
        }

        i <- 0L
        while( i < maxIter){
            i <- i + 1L

            # min_h ||W*H - A||, s.t. H>=0, for given A and W.
            res <- try(.fcnnls(W, A), silent = TRUE)
            if(inherits(res, "try-error")) {
                nrestart <- nrestart+1
                if ( nrestart >= 10 ){
                    warning(paste("als_nmf: Factorization failed too many times",
                                  "[Computation stopped after the 9th",
                                  "restart]"))
                    break
                }
                restart()
                next
            }
            H <- res[[1]]

            if ( any(rowSums(H)==0) ){
                nrestart <- nrestart+1;
                if ( nrestart >= 10 ){
                    warning(paste("als_nmf: Factorization failed too many times",
                                  "[Computation stopped after the 9th",
                                  "restart]"))
                    break
                }
                restart()
                next
            }

            # min_w ||H' * W' - A'||, s.t. W>=0, for given A and H.
            res = try(.fcnnls(t(H), t(A)),
                      silent = TRUE)
            if(inherits(res, "try-error")) {
                nrestart <- nrestart+1;
                if ( nrestart >= 10 ){
                    warning(paste("als_nmf: Factorization failed too many times",
                                  "[Computation stopped after the 9th",
                                  "restart]"))
                    break;
                }
                restart()
                next
            }

            Wt = res[[1]]
            W <- t(Wt);

            #### Convergence test adapted from:####
            ###% M.W. Berry et al. (2007), "Algorithms and Applications for
            ###% Approximate Nonnegative Matrix Factorization," Computational
            ###% Statistics and Data Analysis, vol. 52, no. 1, pp. 155-173.
            # According to the Matlab nnmf, uses RMS instead of Frobenius Norm.
            # Also stops if W and H have converged
            residNorm <- sqrt(mean((A - W %*% H) ^ 2))
            dW <- W - Wold
            dW <- sqrt(mean(dW ^ 2))
            dH <- H - Hold
            dH <- sqrt(mean(dH ^ 2))

            if(abs(residNormOld - residNorm) <= eps_conv ||
               (dW <= eps_conv && dH <= eps_conv)) {
                if(verbose) {
                    cat("Track:\tIter\tNorm\t\trms(dW)\t\trms(dH)\n")
                    cat(sprintf("\t%d\t%f\t%f\t%f\n",
                                i,residNorm, dW, dH))
                    message("Converged!")
                }
                break
            }
            residNormOld <- residNorm
            # end adapted

            Hold <- H
            Wold <- W

            # every 10 iterations
            if (i %% 10==0){
                if ( verbose ){ # prints number of changing elements
                    cat("Track:\tIter\tNorm\t\trms(dW)\t\trms(dH)\n")
                    cat(sprintf("\t%d\t%f\t%f\t%f\n",
                                i,residNorm, dW, dH))
                }
            }
        }
        return(list(obj = residNorm, W = W, H = H))
    })

    solution.best <- which.min(vapply(solutions, function(x) x$obj, numeric(1)))

    W <- solutions[[solution.best]]$W
    H <- solutions[[solution.best]]$H

    res <- new("genericFactorization", W= W, H = H)
    res <- new("genericFit", fit = res, method = "als-nmf")

    return(res)
}

nipals_pca_nocatch <- function(A, k, duplicable = TRUE, ...) {
    if(duplicable) {
        oldSeed <- duplicable("biclus") # do not modify the R global environment
        on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    }

    args <- list(...)
    args <- c(scale = args$scale, center = args$center, maxiter = args$maxiter,
              tol = args$tol,
              startcol = args$startcol, fitted = args$fitted,
              force.na = args$force.na, gramschmidt = args$gramschmidt,
              verbose = args$verbose)
    np <- tryCatch({
        do.call(nipals::nipals, c(list(x = A, ncomp = k,
                                       tol = 1e-6), args))
    },
    error = function(e) {
        if(grepl(pattern = paste0("replacement has length zero"), x = e)) {
            stop(paste("NIPALS was not able to run. Try running clean() on",
                       "your input matrix or BiclusterExperiment object."))
        } else {
            stop(e)
        }
    })

    new("genericFit", fit = new("genericFactorization",
                                W = np$scores,
                                H = t(np$loadings)),
        method = "nipals-pca")
}

#' Principal component dimensionality reduction using NIPALS
#'
#' Factorizes matrix \code{A} as the product of score and loading matrices
#' respectively truncated to \code{k} rows and \code{k} columns. Uses the
#' Nonlinear Iterative Partial Least Squares algorithm to compute principal
#' components in the presence of missing matrix elements.
#'
#' If NIPALS fails, this function will recursively call itself with decreasing
#' values of \code{cleanParam} until NIPALS succeeds.
#'
#' @param A the matrix to factorize
#' @param k the number of factors to compute
#' @param cleanParam passed to \code{\link{clean}()}
#' @param duplicable fix the random seed internally
#' @param verbose report recursive calls and all values of \code{cleanParam}
#' @param ... Additional parameters will be passed to
#'   \code{\link[nipals]{nipals}}.
#'
#' @return a list containing
#'   \describe{
#'   \item{m}{the data matrix after any cleaning}
#'   \item{genericFit}{a \code{\link{genericFit-class}} object}
#'   \item{indexRemaining}{a list of the row and column indexes remaining after
#'   cleaning}
#'   }
#' @export
nipals_pca <- function(A, k, cleanParam = 0,
                       duplicable = FALSE, verbose = TRUE, ...) {
    if(duplicable) {
        oldSeed <- duplicable("biclus") # do not modify the R global environment
        on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    }
    cleanRes <- clean(A, cleanParam, dimsRemain = TRUE)
    mClean <- cleanRes$obj
    indexRem <- cleanRes$dimsRemain

    tryCatch({
        return(list(m = mClean, genericFit = nipals_pca_nocatch(mClean, k, ...),
                    indexRemaining = indexRem))
    }, error = function(e) {
        if(grepl(pattern = paste0("replacement has length zero"), x = e)) {
            cleanParam <- cleanParam + min(1, log10(2.5 - cleanParam))
            if(verbose) {
                message(paste("Too many NA in the data. Cleaning with",
                              "cleanParam at", cleanParam))
            }
            # pass the original m so indexRemaining is valid for the user's matrix
            return(nipals_pca_nocatch(A, k, cleanParam, duplicable))
        } else { stop(e) }
    })
}

# Returns a biclust::Biclust-class object
plaid <- function(A, k, duplicable = TRUE, verbose = TRUE, ...) {
    if(duplicable) {
        oldSeed <- duplicable("biclus") # do not modify the R global environment
        on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    }

    args <- list(...)
    row.release <- list(...)$row.release
    col.release <- list(...)$col.release
    args <- c(cluster = args$cluster, shuffle = args$shuffle,
              fit.model = args$fit.model,
              background = args$background,
              background.layer = args$background.layer,
              background.df = args$background.df,
              backfit = args$backfit, iter.startup = args$iter.startup,
              iter.layer = args$iter.layer, verbose = args$verbose)
    number <- 0
    if(is.null(row.release)) { row.release <- 1 }
    if(is.null(col.release)) { col.release <- 1 }
    best <- NULL
    while(number < k && (row.release > 0 || col.release > 0)) {
        dummy <- capture.output({
            bc <- do.call(biclust::biclust,
                          c(list(x = A, method = biclust::BCPlaid(),
                                 row.release = row.release,
                                 col.release = col.release,
                                 max.layers = k), args))
        })
        if(bc@Number > number) {
            number <- bc@Number
            best <- bc
        }
        if(verbose) {
            cat(paste("method =", class(bc@Parameters$Call$method), "\n"))
            cat(paste("max.layers =", bc@Parameters$Call$max.layers, "\n"))
            cat(paste("row.release =", bc@Parameters$Call$row.release, "\n"))
            cat(paste("col.release =", bc@Parameters$Call$col.release, "\n"))
            cat(paste("Biclusters:", bc@Number, "\n"))
        }

        row.release <- max(0, row.release - 0.1)
        col.release <- max(0, col.release - 0.1)
        # release decrements from 0.7 to 0.1
    }

    if(k > number) {
        k <- number
        warning(paste("Plaid could only find", k, "biclusters"))
    }

    return(best)
}

snmf <- function(A, k, verbose = TRUE, duplicable = TRUE, ...) {
    if(duplicable) {
        oldSeed <- duplicable("biclus") # do not modify the R global environment
        on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    }

    eta <- list(...)$eta
    beta <- list(...)$beta
    args <- list(...)
    args <- c(maxIter = args$maxIter, bi_conv = args$bi_conv,
              eps_conv = args$eps_conv, .options = args$.options)
    res <- NULL

    number <- 0

    if(is.null(beta)) { beta <- mean(A) }
    if(is.null(eta)) {eta <- mean(A) }
    if(eta < 0 || beta < .Machine$double.eps) {
        stop("eta must be >= 0 and beta must be > 0")
    }

    while(eta >= 0 && beta >= .Machine$double.eps && !inherits(res, "NMFfit")) {
        res <- tryCatch({
            # /r causes columns (sample membership) to be sparse.
            # Presumably, samples are more likely than features to
            # correspond 1:1 with a bicluster
            res <- suppressMessages(
                do.call(NMF::nmf, c(list(x = A, rank = k,
                                         method = "snmf/r",
                                         beta = beta, eta = eta,
                                         rng = .Random.seed), args))
            )

            if(verbose) {
                cat(paste("method:", res@method, "\n"))
                cat(paste("parameters:\nbeta:", res@parameters[[1]], "\neta:",
                          res@parameters[[2]], "\n"))
            }
            res
        },
        warning = function(w) {
            if (any(suppressWarnings(
                grepl(
                    "too big 'beta' value",
                    w$message,
                    fixed = TRUE
                )
            ))) {
                beta <<- beta / 2

                if(verbose) {
                    cat(paste("NMF::nmf:", w$message, "\n"))
                    cat(paste("mfBiclust::snmf: Restarting with beta = ", beta,
                              "\n"))
                }
            } else {
                warning(w)
            }
            NULL
        },
        error = function(e) {
            if(any(suppressWarnings(grepl("system is computationally singular",
                         e[[1]], fixed = TRUE)))) {
                # While this error is thrown, try decreaseing eta
                if(eta > .Machine$double.eps) {
                    eta <<- eta / 2
                } else {
                    eta <<- 0
                }
                if(verbose) {
                    cat(paste("NMF::nmf:", e[[1]], "\n"))
                    cat(paste("mfBiclust::snmf: Restarting with eta = ", eta,
                              "\n"))
                }
            } else if(any(suppressWarnings(grepl("system is exactly singular",
                                                 e[[1]], fixed = TRUE)))) {
                # When this error thrown, just abort.
                if(verbose) {
                    cat(paste("NMF::nmf:", e[[1]], "\n"))
                }
                stop(message = "Aborting snmf because dataset is singular\n",
                     call. = TRUE)
                # FIXME causes error in biclusterGUI
                # res <- try({
                #     w <- cbind(A[ , 1], matrix(0, nrow = nrow(A), ncol = (k - 1)))
                #     h <- MASS::ginv(w) %*% A
                #     res <- new("genericFactorization", W = w, H = h)
                #     new("genericFit", fit = res, method = "snmf")
                # }, silent = TRUE)
            }
            if(inherits(res, "genericFit"))
                return(res)
            else
                return(NULL)
        })
    }
    if(inherits(res, "NMFfit") || inherits(res, "genericFit")) {
        res@method <- "snmf"
        return(res)
    } else {
        res <- new("genericFactorization", W = matrix(data = numeric()),
                   H = matrix(numeric()))
        res <- new("genericFit", fit = res, method = "snmf")
        return(res)
    }
}

# Returns a biclust::Biclust-class object
spectral <- function(A, k, minSize = NULL, reps = 1, duplicable = TRUE,
                     verbose = TRUE, ...) {
    # spectral may find over k biclusters, but only k will be returned

    # minSize can be used to force biclusters to be a certain fraction of the
    # smaller matrix dimension
    if(duplicable) {
        oldSeed <- duplicable("biclus") # do not modify the R global environment
        on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
    }

    minDim <- min(nrow(A), ncol(A))
    minx <- if(is.null(minSize)) 2 else floor(minDim * minSize)
    if(nrow(A) < 6 || ncol(A) < 6) {
        stop("For Spectral the minimum size of m is 6x6")
    }

    # sometimes withinVar is explicitly named by the caller
    withinVar <- list(...)$withinVar

    # JNL the strategy is to find the lowest value of withinVar that yields
    # enough biclusters.

    # withinVar = minDim is an arbitrary starting point used in some comparison
    # studies
    if(is.null(withinVar)) {
        withinVar <- minDim
    }

    numberOfEigenvalues <- list(...)$numberOfEigenvalues
    if(is.null(numberOfEigenvalues)) {
        # the number of eigenvalues to consider. This can limit number of
        # biclusters found, so increase from the default if necessary.
        numberOfEigenvalues <- 3

        while(!(floor(min(100, nrow(A) / minx) - 1) * # the max number of row clusters
                floor(ncol(A) / minx - 1)) * # the max number of col clusters
              sum(seq_len(numberOfEigenvalues)) >= k) {
            # internally, every row cluster i in e is combined with every col
            # cluster >= i so the number of possible biclusters is the summation
            # of numbers up to e
            numberOfEigenvalues <- numberOfEigenvalues + 1
        }
    }

    number <- 0 # save the biclustering solution with the most clusters
    best <- NULL
    while(number < k) {
        bc <- do.call(
            biclust::biclust, list(x = A, method = biclust::BCSpectral(),
                                   normalization = "log", withinVar= withinVar,
                                   minr = minx, minc = minx,
                                   numberOfEigenvalues = numberOfEigenvalues))
        if(bc@Number > number) {
            number <- bc@Number
            best <- bc
        }
        withinVar <- withinVar + minDim

        if(verbose) {
            cat(paste("method =", class(bc@Parameters$Call$method), "\n"))
            cat(paste("normalization =", bc@Parameters$Call$normalization, "\n"))
            cat(paste("withinVar =", bc@Parameters$Call$withinVar, "\n"))
            cat(paste("minr =", bc@Parameters$Call$minr, "\n"))
            cat(paste("minc =", bc@Parameters$Call$minc, "\n"))
            cat(paste("numberOfEigenvalues =", bc@Parameters$Call$numberOfEigenvalues, "\n"))
            cat(paste("Biclusters:", bc@Number, "\n"))
        }
    }
    if(k > number) {
        k <- number
        warning(paste("Spectral could only find", k, "biclusters"))
    }

    return(best)
}

#' Principal component dimensionality reduction using SVD
#'
#' Factorizes matrix \code{A} as the product of score and loading matrices
#' respectively truncated to \code{k} rows and \code{k} columns.
#'
#' @param A the matrix to factorize
#' @param k the number of factors to compute
#' @param duplicable fix the random seed internally
#' @param ... Present for compatibility. \code{als_nmf} has no other
#'   parameters
#'
#' @return a \code{\link{genericFit-class}} object
#' @export
svd_pca <- function(A, k, duplicable = NULL, ...) {
    prcmp <- prcomp(t(A), rank. = k, retx = TRUE, center = FALSE)
    new(
        "genericFit",
        fit = new(
            "genericFactorization",
            H = t(prcmp$x),
            W = prcmp$rotation
        ),
        method = "svd-pca"
    )
}
