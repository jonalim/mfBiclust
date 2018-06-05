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
  
  if (bicluster == "snmf/l" || bicluster == "snmf") {
    # Use NMF package
    tryCatch(bc <- NMF::nmf(m, k, method = bicluster),
             error = function(c) {warning(paste0("Switching to PCA, the ",
                                                 "preferred method for a ",
                                                 "matrix containing negative ",
                                                 "values."))
               prcmp <- prcomp(m, rank. = k, retx = TRUE)
               bc <<- new("genericFit", fit = new("genericFactorization", W = prcmp$x, H = prcmp$rotation), 
                          method = "pca")
               bicluster <<- "pca"
             }
    )
  }
  else if (bicluster == "pca") {
    # use R pca.
    prcmp <- prcomp(m, rank. = k, retx = TRUE)
    bc <- new("genericFit", fit = new("genericFactorization", W = t(prcmp$x), H = prcmp$rotation), 
              method = "pca")
  }
  
  else if (bicluster == "plaid" || bicluster == "bimax") {
    # FIXME use Biclust algorithms
  }
  
  biclustNames <- unlist(sapply(seq_len(ncol(bc@fit@W)), function(x) {
    paste0("Bicluster.", x)
  }
  ))
  colnames(bc@fit@W) <- biclustNames
  colnames(bc@fit@H) <- biclustNames
  
  #### Thresholding ############################################################
  st <- matrix()
  lt <- matrix()
  sta <- ""
  lta <- ""
  if (inherits(bc, "NMFfit") || inherits(bc, "genericFit")) {
    # thresholding needed if matrix-factorization is being performed
    #### Score thresholds ####
    if (is.numeric(scoreThresh)) {
      # User provided threshold(s)
      st <- scoreThresh
      if (inherits(scoreThresh, "matrix")) {
        names <- unlist(lapply(
          seq_len(ncol(scoreThresh)),
          function(x) paste("User no.", bcs, sep = "")
        ))
      }
      else {
        names <- c("User")
      }
    } else {
      # Calculate thresholds using available algorithms
      st <- apply(bc@fit@W, 2, function(x) {
        rescaled <- (x - min(x)) / (max(x) - min(x))
        thresholds <- c(EBImage::otsu(as.matrix(rescaled)))
        thresholds * (max(x) - min(x)) + min(x)
      })
      colNames <- c("otsu")
    }
    # combine threshold values and names into a matrix
    st <- generateThresholdMatrix(st, k, biclustNames, colNames)
    
    # Note that the user wants the otsu threshold to be applied (others can
    # be plotted simultaneously in the GUI)
    if (length(scoreThresh) == 1 && scoreThresh %in% c("otsu")) {
      sta <- scoreThresh
    }
    
    #### Loading thresholds ####
    if (is.numeric(loadingThresh)) {
      # User provided threshold(s)
      lt <- loadingThresh
      if (inherits(loadingThresh, "matrix")) {
        names <- unlist(lapply(
          seq_len(ncol(loadingThresh)),
          function(x) paste("User no.", bcs, sep = "")
        ))
      } else {
        names <- c("User")
      }
    } else {
      # Calculate thresholds using available algorithms
      lt <- apply(bc@fit@H, 2, function(x) {
        rescaled <- (x - min(x)) / (max(x) - min(x))
        thresholds <- c(EBImage::otsu(as.matrix(rescaled)))
        thresholds * (max(x) - min(x)) + min(x)
      })
      names <- c("otsu")
      
      # Note that the user wants the otsu threshold to be applied (others can
      # be plotted simultaneously in the GUI)
      if (length(loadingThresh) == 1 && loadingThresh %in% c("otsu")) {
        lta <- loadingThresh
      }
    }
    
    # combine threshold values and names into a matrix
    lt <- generateThresholdMatrix(lt, k, biclustNames, colNames)
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
# Combines threshold values and names into a matrix
# Warning: does not do input format checking. Names must be same length as the
# expected number of algorithms
generateThresholdMatrix <- function(ts, k, rowNames, colNames) {
  if (inherits(ts, "matrix") && mode(ts) == "numeric" && nrow(ts) == k) {
    # A matrix of numerics, if k x Y for any Y, will be assumed to be a matrix
    # of thresholds, where each row k contains multiple thresholds to plot for
    # bicluster k.
    tMatrix <- ts
    if (!is.null(colnames(ts))) {
      colnames(ts) <- colNames
    }
    if (!is.null(rownames(ts))) {
      rownames(ts) <- rowNames
    }
  } else if (inherits(ts, "numeric") && length(ts) == 1L) {
    # A single numeric will be applied to all clusters
    tMatrix <- matrix(rep(ts, times = k),
                      ncol = 1,
                      dimnames = list(rowNames, colNames)
    )
  } else if (inherits(ts, "numeric") && length(ts) == k) {
    # A vector of numerics, if the same size as k, will be assumed to have
    # a 1:1 relation with k
    tMatrix <- matrix(ts,
                      ncol = 1,
                      dimnames = list(rowNames, colNames)
    )
  } else {
    stop("The format, dimensions, or length of the argument \"ts\"
is incorrect. Please ensure \"ts\" is numeric, and either
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
