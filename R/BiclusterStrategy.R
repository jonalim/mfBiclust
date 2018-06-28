#' @include helperFunctions.R
#' @include containerGenerics.R
NULL

#### CLASS #####################################################################
#threshold algos must be characters pointing to columns of scoreThresh and
#loadingThresh. If NULL, then all functions assume the first column determines
#bicluster membership.
setClass(
  "BiclusterStrategy",
  slots = list(
    factors = "ANY",
    scoreThresh = "matrix",
    loadingThresh = "matrix",
    scoreThreshAlgo = "character",
    loadingThreshAlgo = "character",
    pred = "matrix",
    name = "character"
  )
)

#### CONSTRUCTOR ###############################################################
#' Construct a BiclusterStrategy
#'
#' This class encapsulates bicluster results for one biclustering algorithm, one
#' thresholding algorithms, and one quantity of biclusters. To visualize results
#' in a GUI, wrap a \code{\link{BiclusterStrategy}} in a
#' \code{\link{BiclusterExperiment}}, then call \code{\link{shinyStart()}}.
#'
#' details
#' @section Custom thresholds:
#' When giving custom thresholds, various common use cases are assumed based on
#' the data type: A single numeric will be applied to all clusters. A vector of
#' numerics, if the same size as k, will be assumed to have a 1:1 relation with
#' k. A matrix of numerics, if k x Y for any Y, will be assumed to be a matrix
#' of thresholds, where each row k contains multiple thresholds to plot for
#' bicluster k. The first threshold will be applied to determine bicluster
#' members.
#' 
#' Careful! Use of NIPALS-PCA is not allowed when data is not missing. (NIPALS
#' is very slow.) Note: when Spectral, can provide minR as the minimum fraction of rows/columns that can constitute a cluster.
#'
#' To be added: factorize: sparsenmf, plaid, bimax.
#' To be added: threshold: ita, fcm
#'
#' @param bicluster the biclustering algorithm to use
#' @param scoreThresh the score thresholding algorithm to use. Ignored if
#' bicluster is "plaid" or "bimax"
#' @param loadingThresh the loading thresholding algorithm to use. Ignored if
#' bicluster is "plaid" or "bimax"
BiclusterStrategy <-
  function(m,
           k,
           method = c("als-nmf", "svd-pca", "snmf", "nipals-pca", "plaid", "spectral"),
           scoreThresh = c("otsu"),
           loadingThresh = c("otsu"), ...) {
    method = match.arg(method)
    if (!"matrix" %in% class(m)) {
      warning(paste0(
        "Argument \"m\" must be of type matrix. Attempting to",
        "coerce m to matrix."
      ))
      m <- as.matrix(m)
    }
    
    #### Matrix factorization ###################################################
    bc <- NULL
    
    if(!any(is.na(m))) {
      
      # These three are only valid when m is complete
      if (method == "svd-pca") {
        # use R pca.
        bc <- svd_pca(m, k)
      } else if (method == "als-nmf") {
        if(any(m < 0)) {
          warning(paste("Converting to pseudovalues (x + abs(min(x))) just for",
                        "this BiclusterStrategy because",
                        "negative values are not allowed."))
          m <- m + abs(min(m))
        }
        tryCatch(
          # Use helper function adapted from the NMF package
          bc <- als_nmf(m, k),
          error = function(c) {
            warning("ALS-NMF failed, switching to PCA.")
            bc <<- svd_pca(m, k) # fallback to PCA
            method <<- "svd-pca"
            return(NULL)
          }
        )
      } else if (method == "snmf") {
        # Use NMF package
        tryCatch(
          bc <- snmf(m, k),
          error = function(c) {
            warning("Sparse NMF failed, switching to PCA.")
            bc <<- svd_pca(m, k) # fallback to PCA
            method <<- "svd-pca"
            return(NULL)
          }
        )
      } else if (method == "plaid") {
        bc <- plaid(m, k)
      } else if (method == "spectral") {
        bc <- spectral(m, k, ...)
      }
    } else {
      if (method != "nipals-pca") {
        warning(paste("Switching to the NIPALS-PCA method because the input",
                      "matrix has missing data"))
        method <- "nipals-pca"
      }
      bc <- nipals_pca(m, k)
    }

    k <- ncol(bc@fit@W) # sometimes the biclustering method returns less than
    # k biclusters
    biclustNames <- unlist(sapply(seq_len(k), function(x) {
      paste0("Bicluster.", x)
    }))
    colnames(bc@fit@W) <- biclustNames
    rownames(bc@fit@H) <- biclustNames
    #### Thresholding ############################################################
    st <- matrix()
    lt <- matrix()
    sta <- ""
    lta <- ""
    if (inherits(bc, "NMFfit") || inherits(bc, "genericFit")) {
      # thresholding needed if matrix-factorization is being performed
      #### Score thresholds ####
      st <-
        generateThresholdMatrix(scoreThresh, bc@fit@W, biclustNames)
      lt <-
        generateThresholdMatrix(loadingThresh, t(bc@fit@H), biclustNames)
      
      sta <- colnames(st)[1]
      lta <- colnames(lt)[1]
    } else {
      # leave thresholds NULL? still to be implemented
    }
    
    #### Results #############################################################
    pred <- threshold(bc@fit@W, st[1])
    
    bcs <- new(
      "BiclusterStrategy",
      factors = bc,
      scoreThresh = st,
      loadingThresh = lt,
      scoreThreshAlgo = sta,
      loadingThreshAlgo = lta,
      pred = pred
    )
    bcs@name <- name(bcs)
    bcs
  }

#### METHODS ###################################################################
validBiclusterStrategy <- function(object) {
  msg <- NULL
  factors <- object@factors
  if (!(inherits(factors, "NMFfit") ||
        inherits(factors, "genericFit"))) {
    msg <-
      c(msg,
        paste(
          "The factors slot must be an 'NMFfit' or 'genericFit'",
          "object"
        ))
  } else {
    val <- validObject(factors, test = TRUE)
    if (inherits(val, "character")) {
      msg <- c(msg, val)
    }
  }

  sThresh <- scoreThresh(object)
  if (!(inherits(sThresh, "matrix") &&
        mode(sThresh) == "numeric")) {
    msg <- c(msg, "The scoreThresh slot must be a numeric matrix")
  } else {
    if (nrow(sThresh) != nclust(object)) {
      msg <-
        c(
          msg,
          paste(
            "The scoreThresh matrix must have number of rows",
            "equal to the width of the score matrix"
          )
        )
    }
    if (!setequal(rownames(sThresh), colnames(score(object)))) {
      msg <-
        c(
          msg,
          paste(
            "scoreThresh column names must correspond 1:1 with",
            "score matrix row names. These are bicluster names"
          )
        )
    }
  }
  
  lThresh <- loadingThresh(object)
  if (!(inherits(lThresh, "matrix") &&
        mode(lThresh) == "numeric")) {
    msg <- c(msg, "The loadingThresh slot must be a numeric matrix")
  } else {
    if (nrow(lThresh) != nclust(object)) {
      msg <-
        c(
          msg,
          paste(
            "The loadingThresh matrix must have number of rows",
            "equal to the width of the score matrix"
          )
        )
    }
    if (!setequal(rownames(lThresh), rownames(loading(object)))) {
      msg <-
        c(
          msg,
          paste(
            "loadingThresh column names must correspond 1:1 with",
            "loading matrix row names. These are bicluster names"
          )
        )
    }
  }
  
  if (!inherits(object@scoreThreshAlgo, "character")) {
    msg <- c(msg, "The scoreThreshAlgo slot must be a character string")
  }
  if (!inherits(object@loadingThreshAlgo, "character")) {
    msg <-
      c(msg, "The loadingThreshAlgo slot must be a character string")
  }
  
  predM <- pred(object)
  if (!(inherits(predM, "matrix") && mode(predM) == "logical")) {
    msg <- c(msg, "The prediction matrix must be a logical matrix")
  } else {
    if (!identical(dim(predM), c(nrow(score(object)), nclust(object)))) {
      msg <-
        c(
          msg,
          paste(
            "The prediction matrix must have dimensions M x K",
            "where M = the height of score(object) and",
            "K = the width of score(object)."
          )
        )
    }
    if (!identical(dimnames(predM), dimnames(score(object)))) {
      msg <-
        c(
          msg,
          paste(
            "The row and column names of the prediction matrix",
            "must be identical to the row and column names of the",
            "score matrix"
          )
        )
    }
  }
  if (is.null(msg))
    TRUE
  else
    msg
}
setValidity("BiclusterStrategy", validBiclusterStrategy)

#### method ####
#' @export
setGeneric("method", signature = "bcs", function(bcs) {
  standardGeneric("method")
})
setMethod("method", "BiclusterStrategy", function(bcs) {
  bcs@factors@method
})

#' Name of a BiclusterStrategy
#'
#' Get this BiclusterStrategy's display-friendly name. If it does not have a name, computes a
#' string containing the biclustering algorithm, thresholding algorithms, and
#' number of biclusters in the given BiclusterStrategy.
#'
#' This function may not be used to modify a BiclusterStrategy's name.
setMethod("name", c(bcs = "BiclusterStrategy"), function(bcs) {
  if (length(bcs@name) > 0) {
    bcs@name
  }
  else {
    # Capitalize 
    # als-nmf is not a valid name for a switch case
    bca <- capitalize(method(bcs))
    sta <- capitalize(bcs@scoreThreshAlgo)
    lta <- capitalize(bcs@loadingThreshAlgo)
  }
  name(list(bca = bca, sta = sta, lta = lta, k = nclust(bcs)))
})

setMethod("name", c(bcs = "list"), function(bcs) {
  paste(
      bcs$bca,
      paste(bcs$sta, bcs$lta, sep = "/"),
      bcs$k,
      sep = " | "
    )
  })

#' Names of biclusters in this BiclusterStrategy
#' @export
setMethod("names", "BiclusterStrategy", function(x) {
  colnames(x@factors@fit@W)
})

setMethod("nclust", c(bcs = "BiclusterStrategy"), function(bcs) {
  ncol(bcs@factors@fit@W)
})
setMethod("nclust", c(bcs = "list"), function(bcs) {
  ncol(bcs[[1]]@factors@fit@W)
})

#' Loading matrix
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is k x N.
#'
#' @export
setGeneric("loading", signature = "bcs", function(bcs) {
  standardGeneric("loading")
})
setMethod("loading", "BiclusterStrategy", function(bcs) {
  bcs@factors@fit@H
})

#' Loading thresholds
#'
#' @export
setGeneric("loadingThresh", signature = "bcs", function(bcs) {
  standardGeneric("loadingThresh")
})
setMethod("loadingThresh", "BiclusterStrategy", function(bcs) {
  bcs@loadingThresh
})

setMethod("pred", c(bcs = "BiclusterStrategy"), function(bcs) {
  bcs@pred
})

#' Score matrix
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is M x k.
#''
#' @export
setGeneric("score", signature = "bcs", function(bcs) {
  standardGeneric("score")
})
setMethod("score", "BiclusterStrategy", function(bcs) {
  bcs@factors@fit@W
})

#' Score thresholds
#'
#' @export
setGeneric("scoreThresh", signature = "bcs", function(bcs) {
  standardGeneric("scoreThresh")
})
setMethod("scoreThresh", "BiclusterStrategy", function(bcs) {
  bcs@scoreThresh
})

#### HELPER FUNCTIONS ##########################################################

#' Combine threshold values and names into a matrix
#'
#' Thresholds may contain a vector of the names of desired threshold algorithms.
#' Alternatively, a numeric-mode object may be supplied. A single numeric
#' will be applied to all clusters. A vector of numerics, if the same size as
#' k, will be assumed to have a 1:1 relation with k. A matrix of numerics,
#' if k x Y, will be assumed to be a matrix where each row k contains Y
#' thresholds to plot for bicluster k. This argument's dimensions will be
#' checked for compatibility with the target matrix.
#'
#' A matrix where each column is a series of thresholds. Each column is
#' named by the threshold algorithm, if specified. Each row is named by the
#' provided biclustNames.
#'
#' @param thresholds either a numeric matrix, a numeric vector, or a character
#'   vector
#' @param matrix the target matrix, whose columns will be thresholded
#' @param biclustNames names of the threshold matrix rows
generateThresholdMatrix <-
  function(thresholds, matrix, biclustNames) {
    if (identical(thresholds, "otsu")) {
      # Calculate thresholds using available algorithms
      tMatrix <- as.matrix(apply(matrix, 2, function(x) {
        if (max(x) != min(x)) {
          rescaled <- (x - min(x)) / (max(x) - min(x))
          m <- as.matrix(rescaled)
          thresholds <- EBImage::otsu(m)
          thresholds * (max(x) - min(x)) + min(x)
        } else  {
          x[1]
        }
        # If all values of x are the same, then the threshold is that value
        # itself
      }), ncol = 1)
      colnames(tMatrix) <- thresholds
    } else if (inherits(thresholds, "matrix") &&
               mode(thresholds) == "numeric" &&
               nrow(thresholds) == k) {
      # A matrix of numerics, if k x Y for any Y, will be assumed to be a matrix
      # of thresholds, where each row k contains multiple thresholds to plot for
      # bicluster k.
      colNames <- unlist(lapply(seq_len(ncol(thresholds)),
                                function(x)
                                  paste("User.", x, sep = "")))
      tMatrix <- thresholds
      if (!is.null(colnames(tMatrix))) {
        colnames(tMatrix) <- colNames
      }
    } else if (inherits(thresholds, "numeric") &&
               length(thresholds) == 1L) {
      # A single numeric will be applied to all clusters
      tMatrix <- matrix(
        rep(thresholds, times = k),
        ncol = 1,
        dimnames = list(biclustNames, paste0("Th.", thresholds))
      )
    } else if (inherits(thresholds, "numeric") &&
               length(thresholds) == k) {
      # A vector of numerics, if the same size as k, will be assumed to have
      # a 1:1 relation with k
      coln <- unlist(sapply(thresholds, paste0("Th.", x)))
      tMatrix <- matrix(thresholds,
                        ncol = 1,
                        dimnames = list(biclustNames, coln))
    } else {
      stop(
        "The format, dimensions, or length of the argument \"thresholds\"
        is incorrect. Please ensure \"thresholds\" is numeric, and either
        atomic, a vector of length k, or an matrix with k rows."
      )
    }
    
    if (!is.null(rownames(tMatrix))) {
      rownames(tMatrix) <- biclustNames
    }
    tMatrix
  }
