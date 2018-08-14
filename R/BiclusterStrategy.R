#' @include generics.R
#' @include helperClasses.R
NULL

#### CLASS #####################################################################
#' Class "BiclusterStrategy" for biclusters
#'
#' This class encapsulates bicluster results for one biclustering algorithm, one
#' thresholding algorithm, and one bicluster quantity \eqn{k}. If the
#' biclustering algorithm returns more than biclusters than the value of
#' \code{k} provided to \code{\link{addStrat}()}, the extra biclusters will
#' be stored with a warning. See note on slot \code{k} below.
#' 
#' @slot factors a \code{\link[NMF]{NMFfit-class}} or \code{\link{genericFit}}
#'   containing a pair of matrices \eqn{L_{m,k}} and \eqn{R_{k,n}}.
#' @slot scoreThresh a vector of thresholds of length \eqn{k}
#' @slot loadingThresh a vector of thresholds of length \eqn{k}
#' @slot threshAlgo the name of the thresholding algorithm, or "user"
#' @slot biclust a \code{\link[biclust]{biclust-class}} object containing
#'   logical bicluster membership matrices and other information
#' @slot k the value of parameter \code{k} in \code{\link{addStrat}()}.
#'   \code{k} causes accessor methods to hide any extra biclusters stored
#'   in this object, unless said accessor methods are called with 
#'   \code{allBc = TRUE}
#' @slot name a display friendly character string describing the object
#' @importClassesFrom biclust Biclust
setClass(
  "BiclusterStrategy",
  slots = list(
    factors = "ANY",
    scoreThresh = "numeric",
    loadingThresh = "numeric",
    threshAlgo = "character",
    biclust = "Biclust",
    k = "integer",
    name = "character"
  )
)

#### CONSTRUCTOR ###############################################################
# The function addStrat() is the exported wrapper that calls this constructor.
#
# @param bicluster the biclustering algorithm to use
# @param scoreThresh the score thresholding algorithm to use. Ignored if
# bicluster is "plaid" or "bimax"
# @param loadingThresh the loading thresholding algorithm to use. Ignored if
# bicluster is "plaid" or "bimax"
setGeneric("BiclusterStrategy", signature = c("obj", "k"),
           function(obj, k, method = c("als-nmf", "svd-pca", "snmf",
                                       "nipals-pca", "plaid", "spectral"),
                    threshAlgo = "otsu", scoreThresh = threshAlgo,
                    loadingThresh = threshAlgo, duplicable = TRUE, ...) 
             standardGeneric("BiclusterStrategy"))
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
             warning(paste(c$message))
             warning(paste(method, "failed, switching to PCA."))
             bc <<- svd_pca(obj, k) # fallback to PCA
             method <<- "svd-pca"
             return(NULL)
           }
  )
  # bc may be NMFfit, biclust::Biclust, or genericFit
  # dispatch to appropriate BiclusterStrategy()
  BiclusterStrategy(obj = bc, k, method, threshAlgo, scoreThresh, loadingThresh, 
                    duplicable, ...)
})
setMethod("BiclusterStrategy", c(obj = "Biclust", k = "numeric"), function(
  obj, k, method, threshAlgo, scoreThresh, loadingThresh, duplicable, ...) {
  
  if(obj@Number < k) {
    k <- obj@Number # sometimes the biclustering method returns less than
  }
  if(ncol(obj@RowxNumber) > k) {
    warning(paste("The biclustering method returned more than k biclusters.",
                  "Use accessor functions with allBc = TRUE to reveal all data."))
  }
  
  biclustNames <- unlist(sapply(seq_len(obj@Number), function(x) {
    paste0("Bicluster.", x)
  }))
  colnames(obj@RowxNumber) <- biclustNames
  rownames(obj@NumberxCol) <- biclustNames
  
  # transfer the unthresholded matrices to a genericFit
  factorz <- new("genericFactorization", W = obj@RowxNumber, H = obj@NumberxCol)
  fit <- new("genericFit", fit = factorz, method = method)
  
  #then threshold the Biclust object
  thRes <- thresholdHelper(bc = obj, k = obj@Number, scoreThresh = scoreThresh,
                           loadingThresh, threshAlgo = threshAlgo)
  
  bcs <- new("BiclusterStrategy", factors = fit, scoreThresh = thRes[[2]],
             loadingThresh = thRes[[3]], threshAlgo = thRes[[1]],
             k = as.integer(k),
             biclust = thRes[[4]])
  bcs@name <- name(bcs)
  bcs
})

setMethod("BiclusterStrategy", c(obj = "genericFit", k = "numeric"), function(
  obj, k, method, threshAlgo, scoreThresh, loadingThresh, duplicable, ...) { 
  fit <- obj
  
  # sometimes the biclustering method returns less than k biclusters
  if(ncol(fit@fit@W) < k) { k <- ncol(fit@fit@W) }
  if(ncol(fit@fit@W) > k) {
    warning(paste("The biclustering method returned more than k biclusters.",
                  "Use accessor functions with allBc = TRUE to reveal all data."))
  }
  
  biclustNames <- unlist(sapply(seq_len(ncol(fit@fit@W)), function(x) {
    paste0("Bicluster.", x)
  }))
  colnames(fit@fit@W) <- biclustNames
  rownames(fit@fit@H) <- biclustNames
  bc <- biclust::BiclustResult(list(Call = fit@method), fit@fit@W, fit@fit@H,
                               ncol(fit@fit@W), list())
  
  thRes <- thresholdHelper(bc = bc, k = ncol(fit@fit@W),
                           scoreThresh = scoreThresh,
                           loadingThresh, threshAlgo = threshAlgo)
  
  bcs <- new("BiclusterStrategy", factors = fit, scoreThresh = thRes[[2]],
             loadingThresh = thRes[[3]], threshAlgo = thRes[[1]], 
             k = as.integer(k),
             biclust = thRes[[4]])
  bcs@name <- name(bcs)
  bcs
})

setMethod("BiclusterStrategy", c(obj = "NMFfit", k = "numeric"), function(
  obj, k, method, threshAlgo, scoreThresh, loadingThresh, duplicable, ...) { 
  
  # assume that NMF methods always return k biclusters
  fit <- obj
  k <- ncol(fit@fit@W)
  
  biclustNames <- unlist(sapply(seq_len(k), function(x) {
    paste0("Bicluster.", x)
  }))
  colnames(fit@fit@W) <- biclustNames
  rownames(fit@fit@H) <- biclustNames
  
  bc <- biclust::BiclustResult(list(Call = fit@method), fit@fit@W, fit@fit@H,
                               k, list())
  
  thRes <- thresholdHelper(bc = bc, k = k, scoreThresh = scoreThresh,
                           loadingThresh, threshAlgo = threshAlgo)
  
  bcs <- new("BiclusterStrategy", factors = fit, scoreThresh = thRes[[2]],
             loadingThresh = thRes[[3]], threshAlgo = thRes[[1]], 
             k = as.integer(k),
             biclust = thRes[[4]])
  bcs@name <- name(bcs)
  bcs
})

#### VALIDITY METHOD ####
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
    if (!setequal(names(sThresh), colnames(fuzzyFeatures(object)))) {
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
    if (!setequal(names(lThresh), rownames(fuzzySamples(object)))) {
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
    if (!identical(dim(predS), dim(fuzzySamples(object)))) {
      msg <-
        c(
          msg,
          paste(
            "clusteredSamples must have dimensions identical to",
            "fuzzySamples(object)."
          )
        )
    }
    if (!identical(dimnames(predS), dimnames(fuzzySamples(object)))) {
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
    if (!identical(dim(predF), dim(fuzzyFeatures(object)))) {
      msg <-
        c(
          msg,
          paste(
            "clusteredFeatures must have dimensions identical to",
            "t(fuzzySamples(object))."
          )
        )
    }
    if (!identical(dimnames(predF), dimnames(fuzzyFeatures(object)))) {
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
  
  res <- validObject(biclust(object), test = TRUE)
  if(inherits(res, "character")) { msg <<- c(msg, res) }
  if (is.null(msg))
    TRUE
  else
    msg
}
setValidity("BiclusterStrategy", validBiclusterStrategy)

#### Accessors ####

setMethod("biclust", "BiclusterStrategy", function(bcs) {
  bcs@biclust
})
setMethod("method", "BiclusterStrategy", function(bcs) {
  bcs@factors@method
})

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

#' Name of a BiclusterStrategy
#'
#' Get this BiclusterStrategy's display-friendly name. If it does not have a name, computes a
#' string containing the biclustering algorithm, thresholding algorithm, and
#' number of biclusters in the given BiclusterStrategy.
#'
#' This function may not be used to modify a BiclusterStrategy's name.
#' @export
setMethod("name", c(bcs = "list"), function(bcs) {
  do.call(paste, c(bcs, list(sep = " | ")))
})

#' Names of biclusters in this BiclusterStrategy
#' @export
setMethod("bcNames", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) {return(colnames(bcs@factors@fit@W)) }
  return(colnames(bcs@factors@fit@W)[seq_len(nclust(bcs))])
})

#' @describeIn BiclusterStrategy Number of clusters in a BiclusterStrategy
setMethod("nclust", c(bcs = "BiclusterStrategy"), function(bcs, allBc) {
  if(allBc) { return(ncol(bcs@factors@fit@W)) }
  return(bcs@k)
})

setMethod("fuzzySamples", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) { return(bcs@factors@fit@H) }
  return(bcs@factors@fit@H[seq_len(nclust(bcs)), , drop = FALSE])
})

setMethod("loadingThresh", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) { return(bcs@loadingThresh) }
  return(bcs@loadingThresh[seq_len(nclust(bcs))])
})

#' @describeIn BiclusterStrategy A binary matrix showing sample-bicluster
#' membership
#' @export
setMethod("clusteredSamples", c(bcs = "BiclusterStrategy"), function(
  bcs, allBc) {
  if(allBc) { return(bcs@biclust@NumberxCol) }
  return(bcs@biclust@NumberxCol[seq_len(nclust(bcs)), , drop = FALSE])
})

setMethod("fuzzyFeatures", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) {return(bcs@factors@fit@W) }
  return(bcs@factors@fit@W[, seq_len(nclust(bcs)), drop = FALSE])
})

setMethod("scoreThresh", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) {return(bcs@scoreThresh)}
  return(bcs@scoreThresh[seq_len(nclust(bcs))])
})

#' @describeIn BiclusterStrategy A binary matrix showing feature-bicluster
#' membership
#' @export
setMethod("clusteredFeatures", c(bcs = "BiclusterStrategy"), function(
  bcs, allBc) {
  if(allBc) { return(bcs@biclust@RowxNumber) }
  return(bcs@biclust@RowxNumber[, seq_len(nclust(bcs)), drop = FALSE])
})


#' @describeIn BiclusterStrategy The algorithm used to calculate biclustering
#' thresholds
#' @export
setMethod("threshAlgo", c(bcs = "BiclusterStrategy"), function(bcs) {
  bcs@threshAlgo
})

#### HELPER FUNCTIONS ##########################################################

#' Calculate Otsu thresholds on each column of a matrix
#' 
#' @param matrix the target matrix, whose columns will be thresholded
otsuHelper <- function(matrix) {
  # Calculate thresholds using available algorithms
  thresholds <- apply(matrix, 2, function(x) {
    if (max(x) != min(x)) {
      rescaled <- (x - min(x)) / (max(x) - min(x))
      m <- as.matrix(rescaled)
      thresholds <- EBImage::otsu(m)
      thresholds * (max(x) - min(x)) + min(x)
    } else  {
      # If all values of x are the same, then the threshold is that value itself
      x[1]
    }
  })
  thresholds
}

biclust2genericFit <- function(biclust) {
  if(!inherits(biclust, "Biclust")) {
    stop(paste("biclust2genericFit must be called on a \"Biclust\" class",
               "object"))
  }
  biclusters <- biclust::biclusternumber(biclust)
  scoreLoading <- if(biclusters > 0) { 
    biclusterNumber2scorefuzzySamples(biclusters, A, k) 
  } else { list(matrix(rep(NA, nrow(A)), ncol = 1), 
                matrix(rep(NA, ncol(A)), nrow = 1)) 
  }
  
  new("genericFit", fit = new("genericFactorization", W = scoreLoading[[1]], 
                              H = scoreLoading[[2]]),
      method = object@Parameters$Call)
}

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
      scoreThresh <- otsuHelper(bc@RowxNumber)
      loadingThresh <- otsuHelper(t(bc@NumberxCol))
    } else {
      stop(paste("Currently \"otsu\" is the only implemented thresholding",
                 "algorithm"))
    }
    #### Results #############################################################
    bc@RowxNumber <- threshold(m = bc@RowxNumber, th = scoreThresh, MARGIN = 2)
    bc@NumberxCol <- threshold(m = bc@NumberxCol, th = loadingThresh,
                               MARGIN = 1)
  } else { 
    bc@NumberxCol <- matrix(dimnames = dimnames(bc@NumberxCol))
    bc@RowxNumber <- matrix(dimnames = dimnames(bc@RowxNumber))
    warning(paste("Biclustering did not find valid results. No samples or",
                  "features are biclustered."))
  }
  return(list(threshAlgo, scoreThresh, loadingThresh, bc))
}


#### threshold ####
#' Apply threshold to a score or loading matrix
#'
#' Returns a binary matrix of the same size as \code{m} where all elements over
#' the threshold are 1.
#'
#' If th is a vector, the first element of th will be used as threshold for the
#' first col/row in m, etc.
#' 
#' @export
setMethod("threshold", c(m = "matrix", th = "numeric"), function(m, th,
                                                                 MARGIN) {
  # Get all values further from 0 than the provided threshold
  if(MARGIN == 1) {
    if(length(th) != nrow(m)) { stop("Length of th must equal nrow(m).") }
    mat <- do.call(rbind, lapply(seq_len(nrow(m)), function(row) {
      if(th[row] < 0) compare <- `<` else compare <- `>` 
      compare(m[row, ], th[row])
    }))
  } else {
    if(length(th) != ncol(m)) { stop("Length of th must equal ncol(m)") }
    mat <- do.call(cbind, lapply(seq_len(ncol(m)), function(col) {
      if(th[col] < 0) compare <- `<` else compare <- `>` 
      compare(m[, col], th[col])
    }))
  }
  colnames(mat) <- colnames(m)
  rownames(mat) <- rownames(m)
  return(mat)
}
)
