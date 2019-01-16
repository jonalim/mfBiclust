#### addStrat ####
#' Add a BiclusterStrategy to a BiclusterExperiment
#'
#' Returns a BiclusterExperiment identical to \code{bce} with an additional
#' BiclusterStrategy accessible using \code{strategies()} or \code{getStrat()}.
#'
#' Argument \code{method} is used to compute a number of biclusters. Matrix
#' factorization methods will store matrix factors in the \code{factors} slot of
#' the \code{BiclusterStrategy}. The matrix factors will be thresholded to yield the
#' binary matrices in the \code{clusteredSamples} and \code{clusteredFeatures}
#' slots of the \code{BiclusterStrategy}.
#'
#' @section Potential side effects: Due to requirements of various biclustering
#'   methods, this function may override user parameters, with warning. Also, if
#'   any elements of \code{abund(BiclusterExperiment)} are missing, the row and
#'   column containing those elements may be removed, with warning.
#'
#' @param bce A \code{\link{BiclusterExperiment-class}} object to analyze
#' @param k The number of biclusters to find
#' @param method The biclustering algorithm
#' @param duplicable Makes biclustering deterministic
#' @param silent Suppresses warnings and messages
#'
#' @return A copy of \code{bce} with a new BiclusterStrategy added to
#'  \code{BiclusterExperiment@@strategies}
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2)
#'
#' @export
setGeneric("addStrat", signature = c("bce", "k"),
           function(bce, k, method = c("als-nmf", "svd-pca", "snmf",
                                       "nipals-pca", "plaid", "spectral"),
                    duplicable = TRUE, silent = FALSE, ...) {
               standardGeneric("addStrat")
})
#' @describeIn addStrat Default method
setMethod(
    "addStrat", c(bce = "BiclusterExperiment", k = "numeric"),
    function(bce, k, method = c("als-nmf", "svd-pca", "snmf", "nipals-pca",
                                "plaid", "spectral"),
             duplicable, silent, ...) {
        # Validate parameters
        # k must be whole number, smaller than both dimensions of m
        m <- as.matrix(bce)
        method <- match.arg(method)
        if(length(method) > 1) {
            stop("Argument \"method\" must be a single string")
        }
        # do this so that the recursive calls with a NULL method are also
        # silent
        method.orig <- method

        k <- validateKM(k, m, method)
        # Special code for NIPALS or missing data
        if (method == "nipals-pca" || any(is.na(m))) {
            if(method != "nipals-pca") {
                warning(paste("Since some data is NA, the NIPALS-PCA",
                              "algorithm must be used."))
            }
            nipals.res <- nipals_pca(A = m, cleanParam = 0,
                                     k = k, center = FALSE,
                                     duplicable = duplicable)
            bcs <- BiclusterStrategy(obj = nipals.res$genericFit, k = k,
                                     method = "nipals-pca")
            oldDims <- dim(bce)
            bce <- bce[unlist(nipals.res$indexRemaining[[1]]),
                       unlist(nipals.res$indexRemaining[[2]])]
            if(!identical(oldDims, dim(bce))) {
                warning(paste("Some samples or features with too much missing",
                              "data were removed. sampleNames(bce) and",
                              "featureNames(bce) can be called to see the",
                              "remaining samples and features."))
            }
        } else {
            #### DEFUALT CALL ####
            bcs <- BiclusterStrategy(obj = m, k = k, method = method,
                                     duplicable = duplicable, verbose = !silent,
                                     ...)
        }

        name <- name(bcs)
        strategies(bce)[[name]] <- bcs
        if(validObject(bce)) {
            message(paste("Added BiclusterStrategy named", name))
            return(bce)
        }
    })
