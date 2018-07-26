#' @include helperFunctions.R
#' @include containerGenerics.R
#' @include helperClasses.R
NULL

#### CLASS #####################################################################
#threshold algos must be characters pointing to columns of scoreThresh and
#loadingThresh. If NULL, then all functions assume the first column determines
#bicluster membership.
setClass(
  "BiclusterStrategy",
  slots = list(
    factors = "ANY",
    scoreThresh = "numeric",
    loadingThresh = "numeric",
    threshAlgo = "character",
    clusteredFeatures = "matrix", # m x k
    clusteredSamples = "matrix", # n x k
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
setGeneric("BiclusterStrategy", signature = c("obj", "k"),
           function(obj, k, method = c("als-nmf", "svd-pca", "snmf",
                                       "nipals-pca", "plaid", "spectral"),
                    threshAlgo = "otsu", scoreThresh = threshAlgo,
                    loadingThresh = threshAlgo, duplicable = TRUE, ...) {
             standardGeneric("BiclusterStrategy")
           })
setMethod("BiclusterStrategy", c(obj = "matrix", k = "numeric"), function(
  obj, k, method, threshAlgo, scoreThresh, loadingThresh, duplicable, ...) {
  method <- match.arg(method)
  
  #### Matrix factorization ###################################################
  bc <- NULL
  
  if(any(is.na(obj)) || method == "nipals-pca") {
    if (method != "nipals-pca") {
      warning(paste("Switching to the NIPALS-PCA method because the input",
                    "matrix has missing data"))
      method <- "nipals-pca"
    }
    # If still too many NAs, an error will be thrown back to addStrat
    bc <- nipals_pca(biclustArgs)
  }
  
  # Simply shift all values to be >= 0, if necessary for the algorithm
  if(method == "als-nmf" || method == "snmf") obj <- pseudovalues(obj)
  
  # Same arguments regardless of algorithm
  biclustArgs <- c(list(A = obj, k = k, duplicable = duplicable), list(...))
  
  # Function names here have underscores instead of hyphens
  hyphen <- regexpr(pattern = "-", text = method)[[1]]
  if(hyphen > 0) substr(method, start = hyphen, stop = hyphen) <- "_"
  
  tryCatch(bc <- do.call(get(method), biclustArgs), # Bicluster
           error = function(c) {
             warning(paste(method, "failed, switching to PCA."))
             bc <<- svd_pca(obj, k) # fallback to PCA
             method <<- "svd-pca"
             return(NULL)
           }
  )
  
  k <- ncol(bc@fit@W) # sometimes the biclustering method returns less than
  # k biclusters
  biclustNames <- unlist(sapply(seq_len(k), function(x) {
    paste0("Bicluster.", x)
  }))
  colnames(bc@fit@W) <- biclustNames
  rownames(bc@fit@H) <- biclustNames
  
  thRes <- thresholdHelper(bc = bc, k = k, scoreThresh = scoreThresh,
                           loadingThresh, threshAlgo = threshAlgo)
  list(threshAlgo, scoreThresh, loadingThresh, clusteredFeatures, clusteredSamples)
  
  bcs <- new("BiclusterStrategy", factors = bc, scoreThresh = thRes[[2]],
             loadingThresh = thRes[[3]], threshAlgo = thRes[[1]],
             clusteredSamples = thRes[[5]], clusteredFeatures = thRes[[4]])
  bcs@name <- name(bcs)
  bcs
})
setMethod("BiclusterStrategy", c(obj = "genericFit", k = "numeric"), function(
  obj, k, method, threshAlgo, scoreThresh, loadingThresh, duplicable, ...) { 
  method <- match.arg(method)
  bc <- obj
  
  k <- ncol(bc@fit@W) # sometimes the biclustering method returns less than
  # k biclusters
  biclustNames <- unlist(sapply(seq_len(k), function(x) {
    paste0("Bicluster.", x)
  }))
  colnames(bc@fit@W) <- biclustNames
  rownames(bc@fit@H) <- biclustNames
  
  thRes <- thresholdHelper(bc = bc, k = k, scoreThresh = scoreThresh,
                           loadingThresh, threshAlgo = threshAlgo)
  list(threshAlgo, scoreThresh, loadingThresh, clusteredFeatures, clusteredSamples)
  
  bcs <- new("BiclusterStrategy", factors = bc, scoreThresh = thRes[[2]],
             loadingThresh = thRes[[3]], threshAlgo = thRes[[1]],
             clusteredSamples = thRes[[5]], clusteredFeatures = thRes[[4]])
  bcs@name <- name(bcs)
  bcs
})

thresholdHelper <- function(bc, k, scoreThresh, loadingThresh, threshAlgo) {
  if(k > 0) {
    #### Thresholding ############################################################
    if(inherits(scoreThresh, "numeric") && inherits(loadingThresh, "numeric")) {
      if(length(scoreThresh) != k || length(loadingThresh) != k) {
        stop("Length of \"scoreThresh\" and \"loadingThresh\" must equal \"k\"")
      }
      threshAlgo <- "user"
    } else if(threshAlgo == "otsu") {
      # Use automatic thresholding (default)
      scoreThresh <- otsuHelper(bc@fit@W)
      loadingThresh <- otsuHelper(t(bc@fit@H))
    } else {
      stop(paste("Currently \"otsu\" is the only implemented thresholding",
                 "algorithm"))
    }
    #### Results #############################################################
    clusteredFeatures <- threshold(m = bc@fit@W, th = scoreThresh, MARGIN = 2)
    clusteredSamples <- threshold(m = bc@fit@H, th = loadingThresh,
                                  MARGIN = 1)
  } else { 
    clusteredSamples <- matrix(dimnames = dimnames(bc@fit@W))
    clusteredFeatures <- matrix(dimnames = dimnames(bc@fit@W))
    warning(paste("Biclustering did not find valid results. No samples or",
                  "features are biclustered."))
  }
  return(list(threshAlgo, scoreThresh, loadingThresh, clusteredFeatures, clusteredSamples))
}

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
  if (!inherits(sThresh, "numeric")) {
    msg <- c(msg, "The scoreThresh slot must be a numeric vector")
  } else {
    if (length(sThresh) != nclust(object)) {
      msg <-
        c(
          msg,
          paste(
            "scoreThresh must be as long as the width of the score matrix"
          )
        )
    }
    if (!setequal(names(sThresh), colnames(score(object)))) {
      msg <-
        c(
          msg,
          paste(
            "scoreThresh names must correspond 1:1 with score matrix row",
            "names. These are bicluster names"
          )
        )
    }
  }
  
  lThresh <- loadingThresh(object)
  if (!inherits(lThresh, "numeric")) {
    msg <- c(msg, "The loadingThresh slot must be a numeric vector")
  } else {
    if (length(lThresh) != nclust(object)) {
      msg <-
        c(
          msg,
          paste(
            "loadingThresh must be as long as the width of the score matrix"
          )
        )
    }
    if (!setequal(names(lThresh), rownames(loading(object)))) {
      msg <-
        c(
          msg,
          paste(
            "loadingThresh names must correspond 1:1 with loading matrix row",
            "names. These are bicluster names"
          )
        )
    }
  }
  
  if (!inherits(object@threshAlgo, "character")) {
    msg <- c(msg, "The scoreThreshAlgo slot must be a character string")
  }
  predS <- clusteredSamples(object)
  if (!(inherits(predS, "matrix") && mode(predS) == "logical")) {
    msg <- c(msg, "clusteredSamples must be a logical matrix")
  } else {
    if (!identical(dim(predS), dim(loading(object)))) {
      msg <-
        c(
          msg,
          paste(
            "clusteredSamples must have dimensions identical to",
            "loading(object)."
          )
        )
    }
    if (!identical(dimnames(predS), dimnames(loading(object)))) {
      msg <-
        c(
          msg,
          paste(
            "The row and column names of clusteredSamples",
            "must be identical to the row and column names of the",
            "score matrix"
          )
        )
    }
  }
  
  predF <- clusteredFeatures(object)
  if (!(inherits(predF, "matrix") && mode(predF) == "logical")) {
    msg <- c(msg, "clusteredFeatures must be a logical matrix")
  } else {
    if (!identical(dim(predF), dim(score(object)))) {
      msg <-
        c(
          msg,
          paste(
            "clusteredFeatures must have dimensions identical to",
            "t(loading(object))."
          )
        )
    }
    if (!identical(dimnames(predF), dimnames(score(object)))) {
      msg <-
        c(
          msg,
          paste(
            "The row and column names of clusteredFeatures",
            "must be identical to the column and row names, respectively, of",
            "the score matrix"
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

#### Accessors ####
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
    bca <- capitalize(method(bcs))
    ta <- capitalize(bcs@threshAlgo)
    name(list(bca, ta, nclust(bcs)))
  }
})

setMethod("name", c(bcs = "list"), function(bcs) {
  do.call(paste, c(bcs, list(sep = " | ")))
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

#' Sample-bicluster clustering matrix
#' 
#' A binary matrix showing which samples are members of which biclusters.
#'
#'@export
setGeneric("clusteredSamples", signature = "bcs", function(bcs) {
  standardGeneric("clusteredSamples")
})
setMethod("clusteredSamples", c(bcs = "BiclusterStrategy"), function(bcs) {
  bcs@clusteredSamples
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

#' Feature-bicluster clustering matrix
#' 
#' A binary matrix showing which features are members of which biclusters.
#' 
#' @export
setGeneric("clusteredFeatures", signature = "bcs", function(bcs) {
  standardGeneric("clusteredFeatures")
}) 
setMethod("clusteredFeatures", c(bcs = "BiclusterStrategy"), function(bcs) {
  bcs@clusteredFeatures
})

#' Algorithm used to calculate threshold
#' 
#' Bicluster membership is determined by a binary thresholding function. The
#' value of the threshold can be either user-provided or calculated by the
#' Otsu algorithm.
#' 
#' @export
setGeneric("threshAlgo", signature = "bcs", function(bcs) {
  standardGeneric("threshAlgo")
}) 
setMethod("threshAlgo", c(bcs = "BiclusterStrategy"), function(bcs) {
  bcs@threshAlgo
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
otsuHelper <- function(matrix) {
  # Calculate thresholds using available algorithms
  thresholds <- apply(matrix, 2, function(x) {
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
  })
  thresholds
}
