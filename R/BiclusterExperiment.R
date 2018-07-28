#' @include BiclusterStrategy.R
#' @include helperFunctions.R
#' @include containerGenerics.R
NULL

#### CLASS #####################################################################
#' Class "BiclusterExperiment" for data and biclusters
#'
#' This class encapsulates data for one or more biclustering runs derived from
#' the same abundance data. Objects can be created using the
#' \code{\link{BiclusterExperiment}} constructor.
#'
#' @slot data Object of class \code{\link{matrix}}. The original data.
#' @slot annot Object of class \code{\link{data.frame}}. Annotations provided by the user
#' @slot strategies A \code{\link{list}} of \code{BiclusterStrategy} objects
#' 
#' @example R/examples/addStrat-biclusterGUI.R
#' @importClassesFrom Biobase eSet
setClass("BiclusterExperiment", slots = list(
  strategies = "list"
), contains = "eSet")

setAs("ExpressionSet", "BiclusterExperiment", function(from) {
  ad <- Biobase::exprs(from)
  
  # remove genes with any NA
  # naIndex <- which(rowSums(is.na(ad)) > 0)
  # from <- from[-naIndex, ]
  # ad <- Biobase::exprs(from)
  
  # Add "abund" matrix and remove "exprs" matrix from the assayData object
  from <- Biobase::assayDataElementReplace(from, "abund", ad, validate = FALSE)
  from <- Biobase::assayDataElementReplace(from, "exprs", NULL)
  bce <- new("BiclusterExperiment", assayData = Biobase::assayData(from), 
             phenoData = Biobase::phenoData(from), 
             featureData = Biobase::featureData(from), strategies = list())
  if(validObject(bce, test = FALSE)) bce
})

#### CONSTRUCTOR ###############################################################
#' Construct a container for biclustering runs
#'
#' Constructs a \code{\link{BiclusterExperiment-class}} object holding data,
#' with a slot for future biclustering runs performed on that data. To run
#' biclustering, use \code{\link{addStrat}()}.
#' 
#' @param m the data matrix defining this BiclusterExperiment. Should have
#'   samples as columns and features as rows.
#'
#' @example R/examples/addStrat-biclusterGUI.R
#' @export
setGeneric("BiclusterExperiment", function(m, bcs = list(), phenoData = Biobase::annotatedDataFrameFrom(m, byrow = FALSE), featureData = annotatedDataFrameFrom(m, byrow = TRUE), pp = FALSE, maxNa = 0.5) {
  standardGeneric("BiclusterExperiment")
})

# eSet objects store assayData: rows are features and columns are samples.
setMethod("BiclusterExperiment", c(m = "matrix"), function(m, bcs, phenoData, featureData, pp, maxNa) {
  if(pp) {
    m <- clean(m, maxNa)
  }
  
  if(!inherits(phenoData, "AnnotatedDataFrame")) {
    phenoData <- AnnotatedDataFrame(phenoData)
  }
  if(!inherits(featureData, "AnnotatedDataFrame")) {
    featureData <- AnnotatedDataFrame(featureData)
  }
  
  ad <- Biobase::assayDataNew(storage.mode = "list")
  ad$abund <- m
  
  if(inherits(bcs, "list")) {
    names(bcs) <- lapply(bcs, function(bcs) {name(bcs)})
  } else if(!is.null(bcs)) {
    bcs <- list(bcs)
    names(bcs) <- name(bcs[[1]])
  } else {
    bcs <- list()
  }
  new("BiclusterExperiment", assayData = ad, phenoData = phenoData,
      featureData = featureData, strategies = bcs)
})

#### METHODS ###################################################################

validBiclusterExperiment <- function( object ) {
  msg <- NULL
  if(!inherits(object, "BiclusterExperiment")) {
    msg <- c(msg, paste("Cannot validate a", class(object), 
                        "as BiclusterExperiment"))
  }
  if(!"abund" %in% Biobase::assayDataElementNames(object)) {
    msg <- c(msg, "The assayData slot must contain a matrix named 'abund'")
  }
  
  if(inherits(object@strategies, "list")) {
    if(!all(names(object@strategies) == sapply(object@strategies, name))) {
      msg <- c(msg, paste("List names of object@strategies are not identical",
                          "to BiclusterStrategy names"))
    } else {
      validBcs <- unlist(sapply(names(object), function(bcs) {
        inherits(getStrat(object, bcs), "BiclusterStrategy")
      }))
      if (!all(validBcs)) {
        msg <- c(msg, "All strategies must be BiclusterStrategy objects.")
      } else {
        sapply(names(object), function(bcs) { # Check validity of all strategies
          res <- validObject(getStrat(object, bcs), test = TRUE)
          if(inherits(res, "character")) { msg <<- c(msg, res) }
        })
      }
    }
  } else {
    msg <- c(msg, "The strategies slot must be a 'list' object")
  }
  if(is.null(msg)) TRUE else msg
}
setValidity("BiclusterExperiment", validBiclusterExperiment)

#### addStrat ####
#' Add a BiclusterStrategy to a BiclusterExperiment
#' 
#' Returns a BiclusterExperiment identical to \code{bce} with the addition of a
#' BiclusterStrategy accessible using \code{strategies()} or \code{getStrat()}.
#'
#' The provided \code{method} is used to compute a number of biclusters,
#' sets comprising both samples and features. Matrix factorization methods will
#' store intermediate data in the \code{factors} slot of the BiclusterStrategy.
#' That intermediate data is then thresholded to yield the matrices in the
#' \code{clusteredSamples} and \code{clusteredFeatures} slots of the
#' BiclusterStrategy.
#'
#' @section Potential side effects:
#' Due to requirements of various biclustering methods, this function may with
#' warning override user parameters. Also, if any elements of 
#' \code{abund(BiclusterExperiment)} are missing, the row and column
#' containing those elements may be removed with warning.
#' 
#' @example R/examples/addStrat-biclusterGUI.R
#' @export
setGeneric("addStrat", signature = c("bce", "k"), function(bce, k, 
                                                           method = c("als-nmf", "svd-pca", "snmf",
                                                                      "nipals-pca", "plaid", "spectral"),
                                                           duplicable = TRUE, silent = FALSE, 
                                                           ...) {
  standardGeneric("addStrat")
})
setMethod("addStrat", c(bce = "BiclusterExperiment", k = "numeric"), 
          function(bce, k, method = c("als-nmf", "svd-pca", "snmf",
                                      "nipals-pca", "plaid", "spectral"), 
                   duplicable, silent, ...) {
            # Validate parameters
            # k must be whole number, smaller than both dimensions of m
            m <- as.matrix(bce)
            method <- match.arg(method)
            if(length(method) > 1) {
              stop("Argument \"method\" must be a single string")
            }
            # do this so that the recursive calls with a NULL method are also silent
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
                                       duplicable = duplicable, ...)
            }
            
            name <- name(bcs)
            strategies(bce)[[name]] <- bcs
            if(validObject(bce)) {
              message(paste("Added BiclusterStrategy named", name))
              return(bce)
            }
          })

setGeneric("as.matrix")
#' @describeIn BiclusterExperiment Get abundance values in a BiclusterExperiment
#' @export
setMethod("as.matrix", "BiclusterExperiment", function(x) {
  Biobase::assayDataElement(x, "abund")
})

setMethod("clean", c(object = "BiclusterExperiment"), function(object,
                                                               cleanParam) {
  if(!(cleanParam <= 1 && cleanParam >= 0)) {
    stop("Arg \"cleanParam\" must be in the range of 0 to 1.")
  }
  results <- clean(as.matrix(object), maxNa, TRUE)
  # [[2]] contains a vector of indexes of the remaining columns
  # [[1]] contains the cleaned matrix itself
  bce <- object[results[[2]][[2]], results[[2]][[1]]]
  strategies(bce) <- bce@strategies
  
  if(validObject(bce, test = FALSE)) return(bce)
})

setGeneric("getStrat", signature = "bce", function(bce, id) {standardGeneric("getStrat")})
#' @describeIn BiclusterExperiment Get a BiclusterStrategy contained by a
#'   BiclusterExperiment, by providing either name or integer index
#' @export
setMethod("getStrat", c(bce = "BiclusterExperiment"), function(bce, id) {
  strategies(bce)[[id]] 
})

#' @describeIn BiclusterExperiment Character names of BiclusterStrategies in this BiclusterExperiment
#' @export
setMethod("names", "BiclusterExperiment", function(x) names(x@strategies))

setGeneric("wipe", signature = "bce", function(bce) {standardGeneric("wipe")})
#' @describeIn BiclusterExperiment Return this BiclusterExperiment with all BiclusterStrategy objects removed
#' @export
setMethod("wipe", c(bce = "BiclusterExperiment"), function(bce) {
  strategies(bce) <- list()
  return(bce)
})

setGeneric("wipeExcept", signature = c("bce"),
           function(bce, bcs) {standardGeneric("wipeExcept")})
#' @describeIn BiclusterExperiment Return this BiclusterExperiment with all
#'   BiclusterStrategy objects removed except \code{bcs}. Argument \code{bcs}
#'   can be passed as a name, integer index, or the BiclusterStrategy itself.
#' @export
setMethod("wipeExcept", c(bce = "BiclusterExperiment"), function(bce, bcs) {
  if(inherits(bcs, "character") || inherits(bcs, "numeric")) {
    bcs <- getStrat(bce, bcs)
  }
  if(any(names(bce) == name(bcs))) {
    bce@strategies <- list()
    bce@strategies[[name(bcs)]] <- bcs
    return(bce)
  } else {
    stop(paste("The given BiclusterStrategy is not contained by the",
               "given BiclusterExperiment. Please create the BiclusterStrategy",
               "first."))
  }
}
)

setGeneric("strategies", signature = "bce", function(bce) {
  standardGeneric("strategies")
})
#' @describeIn BiclusterExperiment Get/set a list of the BiclusterStrategy objects
#'   contained by this BiclusterExperiment.
#' @export
setMethod("strategies", c(bce = "BiclusterExperiment"), function(bce) {
  bce@strategies
})

setGeneric("strategies<-", function(object, value) standardGeneric("strategies<-"))
setReplaceMethod("strategies", signature(object = "BiclusterExperiment",
                                         value = "list"),
                 function(object, value) { 
                   object@strategies <- value 
                   if(validObject(object, test = FALSE)) {
                     object }
                 }
)
