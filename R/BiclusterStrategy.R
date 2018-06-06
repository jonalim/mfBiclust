#' @include helperFunctions.R
#' @include containerGenerics.R
NULL

#### CLASS #####################################################################
#threshold algos must be characters pointing to columns of scoreThresh and
#loadingThresh. If NULL, then all functions assume the first column determines
#bicluster membership.
setClass("BiclusterStrategy", 
         slots = list(factors = "ANY", 
                      biclustAlgo = "character",
                      scoreThresh = "matrix",
                      loadingThresh = "matrix",
                      scoreThreshAlgo = "character",
                      loadingThreshAlgo = "character",
                      pred = "matrix",
                      name = "character"
         )
)

#' Construct a BiclusterStrategy
#'
#' This class encapsulates bicluster results for one biclustering algorithm, one
#' thresholding algorithms, and one quantity of biclusters. To visualize results
#' in a GUI, wrap a \code{\link{BiclusterStrategy}} in a \code{\link{BiclusterExperiment}}, then call
#' \code{\link{shinyStart()}}.
#' 
#' details
#' @section Custom thresholds:
#' When giving custom thresholds, various common use cases are assumed based on the data type:
#' A single numeric will be applied to all clusters.
#' A vector of numerics, if the same size as k, will be assumed to have
#' a 1:1 relation with k.
#' A matrix of numerics, if k x Y for any Y, will be assumed to be a
#' matrix of thresholds, where each row k contains multiple thresholds to plot
#' for bicluster k. The first threshold will be applied to determine bicluster members.
#'
#' To be added: factorize: sparsenmf, plaid, bimax.
#' To be added: threshold: ita, fcm
#'
#' @param bicluster the biclustering algorithm to use
#' @param scoreThresh the score thresholding algorithm to use. Ignored if
#' bicluster is "plaid" or "bimax"
#' @param loadingThresh the loading thresholding algorithm to use. Ignored if
#' bicluster is "plaid" or "bimax"
#'
#' @export
BiclusterStrategy <- function(m, k, bicluster = c("snmf/l", "pca"),
                              scoreThresh = c("otsu"), 
                              loadingThresh = c("otsu")) {

  if (!"matrix" %in% class(m)) {
    warning(paste0("Argument \"m\" must be of type matrix. Attempting to",
                   "coerce m to matrix."))
    m <- as.matrix(m)
  }
  # k must be whole number, smaller than both dimensions of m
  k <- validateKM(k, m)
  
  #### Matrix factorization ################################################### 
  bc <- NULL
  
  if (bicluster == "pca") {
    # use R pca.
    bc <- pcaWrapper(m, k)
  } else if (bicluster == "snmf/l" || bicluster == "snmf" || bicluster == "nmf") {
    # Use NMF package
    tryCatch(bc <- snmfWrapper(m, k),
             error = function(c) {warning(paste0("Sparse NMF failed, switching",
                                                 "to PCA."))
               bc <<- pcaWrapper(m, k) # fallback to PCA
               bicluster <<- "pca"
             }
    )
  }
  
  else if (bicluster == "plaid" || bicluster == "bimax") {
    # FIXME use Biclust algorithms
  }
  
  biclustNames <- unlist(sapply(seq_len(k), function(x) {
    paste0("Bicluster.", x)
  }
  ))
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
    st <- generateThresholdMatrix(scoreThresh, bc@fit@W, biclustNames)
    lt <- generateThresholdMatrix(loadingThresh, t(bc@fit@H), biclustNames)
    
    # Note that the user wants the otsu threshold to be applied (others can
    # be plotted simultaneously in the GUI)
    if (identical(scoreThresh, "otsu")) {
      sta <- scoreThresh
    }
    if(identical(loadingThresh, "otsu")) {
      lta <- loadingThresh
    }
  } else {
    # leave thresholds NULL?
  }
  
  #### Predictions #############################################################
  pred <- threshold(bc@fit@W, st[1])
  colnames(pred) <- lapply(seq_len(ncol(pred)), function(x) {
    paste0("Bicluster.", x)
  }
  )
  
  bcs <- new("BiclusterStrategy",
             factors = bc, biclustAlgo = bicluster, scoreThresh = st,
             loadingThresh = lt, scoreThreshAlgo = sta,
             loadingThreshAlgo = lta, pred = pred)
  bcs@name <- name(bcs)
  bcs
}

#### METHODS ###################################################################

#' Score matrix 
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is M x k.
#''
#' @export
setGeneric("score", signature = "bcs", function(bcs) {standardGeneric("score")})
setMethod("score", "BiclusterStrategy", function(bcs) {
  bcs@factors@fit@W
}
)

#' Loading matrix
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is k x N.
#'
#' @export
setGeneric("loading", signature = "bcs", function(bcs) {standardGeneric("loading")})
setMethod("loading", "BiclusterStrategy", function(bcs) {
  bcs@factors@fit@H
}
)

#' Names of biclusters in this BiclusterStrategy
#' @export
setMethod("names", "BiclusterStrategy", function(x) {
  colnames(x@factors@fit@W)
}
)

#' Name of a BiclusterStrategy 
#' 
#' Get this BiclusterStrategy's name. If it does not have a name, computes a 
#' string containing the biclustering algorithm, thresholding algorithms, and
#' number of biclusters in the given BiclusterStrategy.
#' 
#' This function may not be used to modify a BiclusterStrategy's name.
setMethod("name", c(bcs = "BiclusterStrategy"), function(bcs) {
  if(length(bcs@name) > 0) {bcs@name}
  else {
    paste(bcs@biclustAlgo, 
          paste(bcs@scoreThreshAlgo, bcs@loadingThreshAlgo, sep = "/"),
          ncol(bcs@factors@fit@W), sep = " | ")
  }
}
)

#' Cluster count accessor
#' 
#' Get the number of biclusters created in this BiclusterStrategy.
#' 
#' This function may not be used to modify a BiclusterStrategy's number of biclusters.
setMethod("nclust", c(bcs = "BiclusterStrategy"), function(bcs) {
  ncol(bcs@factors@fit@W)
}
)

setMethod("pred", c(bcs = "BiclusterStrategy"), function(bcs) {
  bcs@pred
}
)

#### HELPER FUNCTIONS ##########################################################

#' Lower beta if nmf throws warning
snmfWrapper <- function(m, k, beta = 0.01) {
  tryCatch(suppressMessages(res <- NMF::nmf(m, k, method = "snmf/l", beta = beta)),
    warning = function(w) {
      if(any(suppressWarnings(grepl("too big 'beta' value", w$message, ignore.case = TRUE, fixed = TRUE)))) {
        beta <<- beta^2
        message(paste0("Decreased beta (sparsity parameter) to ", beta))
        res <<- snmfWrapper(m, k, beta)
      } else {
        warning(w)
      }
    }, 
    error = function(e) {
      stop(e)
    },
    finally = function() {
      res
    }
  )
}

pcaWrapper <- function(m, k) {
  prcmp <- prcomp(m, rank. = k, retx = TRUE)
  new("genericFit", fit = new("genericFactorization", W = prcmp$x, H = t(prcmp$rotation)), 
             method = "pca")
}

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
generateThresholdMatrix <- function(thresholds, matrix, biclustNames) {
# browser() debug here next
  if(identical(thresholds, "otsu")) {
    # Calculate thresholds using available algorithms
        tMatrix <- as.matrix(apply(matrix, 2, function(x) {
          rescaled <- (x - min(x)) / (max(x) - min(x))
          thresholds <- c(EBImage::otsu(as.matrix(rescaled)))
          thresholds * (max(x) - min(x)) + min(x)
        }), ncol = 1, dimnames = list(biclustNames, "otsu"))
  } else if (inherits(thresholds, "matrix") && mode(thresholds) == "numeric" && nrow(thresholds) == k) {
    # A matrix of numerics, if k x Y for any Y, will be assumed to be a matrix
    # of thresholds, where each row k contains multiple thresholds to plot for
    # bicluster k.
    colNames <- unlist(lapply(
          seq_len(ncol(thresholds)),
          function(x) paste("User.", bcs, sep = "")
        ))
    tMatrix <- thresholds
    if (!is.null(colnames(tMatrix))) {
      colnames(tMatrix) <- colNames
    }
    if (!is.null(rownames(tMatrix))) {
      rownames(tMatrix) <- biclustNames
    }
  } else if (inherits(thresholds, "numeric") && length(thresholds) == 1L) {
    # A single numeric will be applied to all clusters
    colNames <- c("User")
    tMatrix <- matrix(rep(thresholds, times = k),
                      ncol = 1,
                      dimnames = list(biclustNames, colNames)
    )
  } else if (inherits(thresholds, "numeric") && length(thresholds) == k) {
    # A vector of numerics, if the same size as k, will be assumed to have
    # a 1:1 relation with k
    tMatrix <- matrix(thresholds,
                      ncol = 1,
                      dimnames = list(biclustNames, colNames)
    )
  } else {
    stop("The format, dimensions, or length of the argument \"thresholds\"
is incorrect. Please ensure \"thresholds\" is numeric, and either
atomic, a vector of length k, or an matrix with k rows.")
  }
  tMatrix
}

# Apply threshold to a score or loading matrix
setGeneric("threshold", signature = "m", function(m, ...) {standardGeneric("threshold")})
#' @export
setMethod("threshold", c(m = "matrix"), function(m, th) {
  m[m >= th] <- TRUE
  m[m < th] <- FALSE
  m
}
)
#' @export
setMethod("threshold", c(m = "Matrix"), function(m, th) {
  m[m >= th] <- TRUE
  m[m < th] <- FALSE
  m
}
)
