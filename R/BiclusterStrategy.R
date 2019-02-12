#' @include generics.R
#' @include helperClasses.R
NULL

#### CLASS #####################################################################
#' Class "BiclusterStrategy" for biclusters
#'
#' This class encapsulates bicluster results created by a call to
#' \code{\link{addStrat}()}.
#'
#' @section The \code{k} slot:
#' The numbers of biclusters requested in the call to \code{addStrat()} is
#' stored in \code{BiclusterStrategy@@k}. A \code{BiclusterStrategy-class}
#' object may store more than \code{k} biclusters, but accessors only return
#' data on the first \code{k} biclusters by default. To obtain data on all
#' stored biclusters, accessors must be called with \code{allBc = TRUE}.
#'
#' @slot factors a \code{\link[NMF]{NMFfit-class}} or
#'   \code{\link{genericFit-class}} containing a pair of matrices \eqn{L_{m,k}}
#'   and \eqn{R_{k,n}}.
#' @slot featureThresh a vector of thresholds of length \eqn{k}
#' @slot sampleThresh a vector of thresholds of length \eqn{k}
#' @slot threshAlgo the name of the thresholding algorithm, or "user"
#' @slot biclust a \code{\link{Biclust-class}} object containing
#'   logical bicluster membership matrices and other information
#' @slot k the number of biclusters to access when \code{allBc = FALSE}.
#' @slot name a display friendly character string describing the object
#'
#' @param bcs A BiclusterSrategy class
#' @param allBc Return data on all stored biclusters

#' @importClassesFrom biclust Biclust
#' @seealso \code{\link{getStrat}}
setClass(
  "BiclusterStrategy",
  slots = list(
    factors = "ANY",
    featureThresh = "numeric",
    sampleThresh = "numeric",
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
# @param featureThresh the score thresholding algorithm to use. Ignored if
# bicluster is "plaid" or "bimax"
# @param sampleThresh the loading thresholding algorithm to use. Ignored if
# bicluster is "plaid" or "bimax"
setGeneric("BiclusterStrategy", signature = c("obj", "k"),
           function(obj, k, method = c("als-nmf", "svd-pca", "snmf",
                                       "nipals-pca", "plaid", "spectral"),
                    threshAlgo = "otsu", featureThresh = threshAlgo,
                    sampleThresh = threshAlgo, ...)
             standardGeneric("BiclusterStrategy"))
setMethod("BiclusterStrategy", c(obj = "matrix", k = "numeric"), function(
  obj, k, method, threshAlgo, featureThresh, sampleThresh, ...) {
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
  biclustArgs <- c(list(A = obj, k = k), list(...))

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
  BiclusterStrategy(obj = bc, k, method, threshAlgo, featureThresh, sampleThresh,
                    ...)
})
setMethod("BiclusterStrategy", c(obj = "Biclust", k = "numeric"), function(
  obj, k, method, threshAlgo, featureThresh, sampleThresh, ...) {

  if(obj@Number < k) {
    k <- obj@Number # sometimes the biclustering method returns less than
  }
  if(ncol(obj@RowxNumber) > k) {
    warning(paste("The biclustering method returned more than k biclusters.",
                  "Use accessor functions with allBc = TRUE to reveal all data."))
  }

  biclustNames <- vapply(seq_len(obj@Number), FUN.VALUE = character(1), FUN = function(x) {
    paste0("Bicluster.", x)
  })
  colnames(obj@RowxNumber) <- biclustNames
  rownames(obj@NumberxCol) <- biclustNames

  # transfer the unthresholded matrices to a genericFit
  factorz <- new("genericFactorization", W = obj@RowxNumber, H = obj@NumberxCol)
  fit <- new("genericFit", fit = factorz, method = method)

  #then threshold the Biclust object
  thRes <- thresholdHelper(bc = obj, k = obj@Number, featureThresh = featureThresh,
                           sampleThresh, threshAlgo = threshAlgo)

  bcs <- new("BiclusterStrategy", factors = fit, featureThresh = thRes[[2]],
             sampleThresh = thRes[[3]], threshAlgo = thRes[[1]],
             k = as.integer(k),
             biclust = thRes[[4]])
  bcs@name <- name(bcs)
  bcs
})

setMethod("BiclusterStrategy", c(obj = "genericFit", k = "numeric"), function(
  obj, k, method, threshAlgo, featureThresh, sampleThresh, ...) {
  fit <- obj

  # sometimes the biclustering method returns less than k biclusters
  if(ncol(fit@fit@W) < k) { k <- ncol(fit@fit@W) }
  if(ncol(fit@fit@W) > k) {
    warning(paste("The biclustering method returned more than k biclusters.",
                  "Use accessor functions with allBc = TRUE to reveal all data."))
  }

  biclustNames <- vapply(seq_len(ncol(fit@fit@W)), FUN.VALUE = character(1), FUN = function(x) {
    paste0("Bicluster.", x)
  })
  colnames(fit@fit@W) <- biclustNames
  rownames(fit@fit@H) <- biclustNames
  bc <- biclust::BiclustResult(list(Call = fit@method), fit@fit@W, fit@fit@H,
                               ncol(fit@fit@W), list())

  thRes <- thresholdHelper(bc = bc, k = ncol(fit@fit@W),
                           featureThresh = featureThresh,
                           sampleThresh, threshAlgo = threshAlgo)

  bcs <- new("BiclusterStrategy", factors = fit, featureThresh = thRes[[2]],
             sampleThresh = thRes[[3]], threshAlgo = thRes[[1]],
             k = as.integer(k),
             biclust = thRes[[4]])
  bcs@name <- name(bcs)
  bcs
})

setMethod("BiclusterStrategy", c(obj = "NMFfit", k = "numeric"), function(
  obj, k, method, threshAlgo, featureThresh, sampleThresh, ...) {

  # assume that NMF methods always return k biclusters
  fit <- obj
  k <- ncol(fit@fit@W)

  biclustNames <- vapply(seq_len(k), FUN.VALUE = character(1), FUN = function(x) {
    paste0("Bicluster.", x)
  })
  colnames(fit@fit@W) <- biclustNames
  rownames(fit@fit@H) <- biclustNames

  bc <- biclust::BiclustResult(list(Call = fit@method), fit@fit@W, fit@fit@H,
                               k, list())

  thRes <- thresholdHelper(bc = bc, k = k, featureThresh = featureThresh,
                           sampleThresh, threshAlgo = threshAlgo)

  bcs <- new("BiclusterStrategy", factors = fit, featureThresh = thRes[[2]],
             sampleThresh = thRes[[3]], threshAlgo = thRes[[1]],
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

  sThresh <- featureThresh(object)
  if (!inherits(sThresh, "numeric")) {
    msg <- c(msg, "The featureThresh slot must be a numeric vector")
  } else {
    if (length(sThresh) != nclust(object)) {
      msg <-
        c(
          msg,
          paste(
            "featureThresh must be as long as the width of the score matrix"
          )
        )
    }
    if (!setequal(names(sThresh), colnames(featureFactor(object)))) {
      msg <-
        c(
          msg,
          paste(
            "featureThresh names must correspond 1:1 with score matrix row",
            "names. These are bicluster names"
          )
        )
    }
  }

  lThresh <- sampleThresh(object)
  if (!inherits(lThresh, "numeric")) {
    msg <- c(msg, "The sampleThresh slot must be a numeric vector")
  } else {
    if (length(lThresh) != nclust(object)) {
      msg <-
        c(
          msg,
          paste(
            "sampleThresh must be as long as the width of the score matrix"
          )
        )
    }
    if (!setequal(names(lThresh), rownames(sampleFactor(object)))) {
      msg <-
        c(
          msg,
          paste(
            "sampleThresh names must correspond 1:1 with loading matrix row",
            "names. These are bicluster names"
          )
        )
    }
  }

  if (!inherits(object@threshAlgo, "character")) {
    msg <- c(msg, "The featureThreshAlgo slot must be a character string")
  }
  predS <- clusteredSamples(object)
  if (!(inherits(predS, "matrix") && mode(predS) == "logical")) {
    msg <- c(msg, "clusteredSamples must be a logical matrix")
  } else {
    if (!identical(dim(predS), dim(sampleFactor(object)))) {
      msg <-
        c(
          msg,
          paste(
            "clusteredSamples must have dimensions identical to",
            "sampleFactor(object)."
          )
        )
    }
    if (!identical(dimnames(predS), dimnames(sampleFactor(object)))) {
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
    if (!identical(dim(predF), dim(featureFactor(object)))) {
      msg <-
        c(
          msg,
          paste(
            "clusteredFeatures must have dimensions identical to",
            "t(sampleFactor(object))."
          )
        )
    }
    if (!identical(dimnames(predF), dimnames(featureFactor(object)))) {
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

# add setters

#' @describeIn BiclusterStrategy The names of biclusters in the
#'   BiclusterStrategy
#' @export
setMethod("bcNames", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) {return(colnames(bcs@factors@fit@W)) }
  return(colnames(bcs@factors@fit@W)[seq_len(nclust(bcs))])
})

setMethod("biclust", "BiclusterStrategy", function(bcs) {
  bcs@biclust
})

#' @describeIn clusteredFeatures method for BiclusterStrategy objects
#'
#' @export
setMethod("clusteredFeatures", c(bcs = "BiclusterStrategy"), function(
  bcs, allBc) {
  if(allBc) { return(bcs@biclust@RowxNumber) }
  return(bcs@biclust@RowxNumber[, seq_len(nclust(bcs)), drop = FALSE])
})

#' @describeIn clusteredSamples method for BiclusterStrategy objects
#'
#' @export
setMethod("clusteredSamples", c(bcs = "BiclusterStrategy"), function(
  bcs, allBc) {
  if(allBc) { return(bcs@biclust@NumberxCol) }
  return(bcs@biclust@NumberxCol[seq_len(nclust(bcs)), , drop = FALSE])
})

#' @describeIn BiclusterStrategy Retrieves \code{bcs@@factors@@fit@@W}, the
#'   matrix factor mapping features to biclusters
setMethod("featureFactor", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) {return(bcs@factors@fit@W) }
  return(bcs@factors@fit@W[, seq_len(nclust(bcs)), drop = FALSE])
})

#' @describeIn BiclusterStrategy Retrieves \code{bcs@@featureThresh}, the threshold(s)
#'   used to determine each bicluster's feature membership
setMethod("featureThresh", "BiclusterStrategy", function(bcs, allBc) {
    if(allBc) {return(bcs@featureThresh)}
    return(bcs@featureThresh[seq_len(nclust(bcs))])
})

#' @describeIn BiclusterStrategy Retrieves \code{bcs@@sampleThresh}, the threshold(s)
#'   used to determine each bicluster's sample membership
setMethod("sampleThresh", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) { return(bcs@sampleThresh) }
  return(bcs@sampleThresh[seq_len(nclust(bcs))])
})

#' @describeIn BiclusterStrategy Retrieves \code{bcs@@factors@@method}, the
#'   biclustering algorithm used to produce \code{bcs}
setMethod("method", "BiclusterStrategy", function(bcs) {
  bcs@factors@method
})

#' @describeIn BiclusterStrategy Gets a display-friendly name of a
#'   BiclusterStrategy that includes its biclustering algorithm, its
#'   thresholding algorithm, and the number of biclusters.
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
#' Helper function (do not call)
#'
#' @rdname name
setMethod("name", c(bcs = "list"), function(bcs) {
  do.call(paste, c(bcs, list(sep = " | ")))
})

#' @describeIn BiclusterStrategy Number of clusters in a BiclusterStrategy
setMethod("nclust", c(bcs = "BiclusterStrategy"), function(bcs, allBc) {
  if(allBc) { return(ncol(bcs@factors@fit@W)) }
  return(bcs@k)
})

#' @describeIn BiclusterStrategy Retrieves \code{bcs@@factors@@fit@@H}, the
#'   matrix factor mapping samples to biclusters
setMethod("sampleFactor", "BiclusterStrategy", function(bcs, allBc) {
  if(allBc) { return(bcs@factors@fit@H) }
  return(bcs@factors@fit@H[seq_len(nclust(bcs)), , drop = FALSE])
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
#'
#' @return numeric vector
#'
#' @export
#' @examples
#' m <- matrix(rnorm(400), 40, 10)
#' otsuHelper(m)
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
  scoreLoading <- if(biclust@Number > 0) {
    list(biclust@RowxNumber, biclust@NumberxCol)
  } else { list(matrix(rep(NA, nrow(biclust@RowxNumber)), ncol = 1),
                matrix(rep(NA, ncol(biclust@NumberxCol), nrow = 1)))
                # FIXME please test if this works
  }

  new("genericFit", fit = new("genericFactorization", W = scoreLoading[[1]],
                              H = scoreLoading[[2]]),
      method = biclust@Parameters$Call)
}

thresholdHelper <- function(bc, k, featureThresh, sampleThresh, threshAlgo) {
  if(k > 0) {
    #### Thresholding ############################################################
    if(inherits(featureThresh, "numeric") && inherits(sampleThresh, "numeric")) {
      if(length(featureThresh) != k || length(sampleThresh) != k) {
        stop("Length of \"featureThresh\" and \"sampleThresh\" must equal \"k\"")
      }
      threshAlgo <- "user"
    } else if(threshAlgo == "otsu") {
      # Use automatic thresholding (default)
      featureThresh <- otsuHelper(bc@RowxNumber)
      sampleThresh <- otsuHelper(t(bc@NumberxCol))
    } else {
      stop(paste("Currently \"otsu\" is the only implemented thresholding",
                 "algorithm"))
    }
    #### Results #############################################################
    bc@RowxNumber <- threshold(m = bc@RowxNumber, th = featureThresh, MARGIN = 2)
    bc@NumberxCol <- threshold(m = bc@NumberxCol, th = sampleThresh,
                               MARGIN = 1)
  } else {
    bc@NumberxCol <- matrix(dimnames = dimnames(bc@NumberxCol))
    bc@RowxNumber <- matrix(dimnames = dimnames(bc@RowxNumber))
    warning(paste("Biclustering did not find valid results. No samples or",
                  "features are biclustered."))
  }
  return(list(threshAlgo, featureThresh, sampleThresh, bc))
}


#### threshold ####
#' @describeIn threshold Apply threshold to a score or loading matrix
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
