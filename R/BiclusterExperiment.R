#' @include BiclusterStrategy.R
#' @include helperFunctions.R
#' @include containerGenerics.R
NULL

#### CLASS #####################################################################
#' Class "BiclusterExperiment" for multiple biclustering results
#'
#' This class encapsulates factorization and thresholding data for one or more
#' biclustering runs. Objects can be created using the
#' \code{\link{BiclusterExperiment}} constructor.
#'
#' @slot data Object of class \code{\link{matrix}}. The original data.
#' @slot annot Object of class \code{\link{data.frame}}. Annotations provided by the user
#' @slot strategies A \code{\link{list}} of \code{BiclusterStrategy} objects
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
#' Perform multiple biclustering runs
#'
#' BiclusterExperiment constructs an object holding data from multiple
#' biclustering runs.
#'
#' This function is useful for constructing one BiclusterExperiment
#' encapsulating results of different pipelines (e.g. comparing NMF with PCA).
#' For comparing results from the same pipeline with differing values of
#' \code{k}, \code{\link{mfbc()}} is recommended.
#' 
#' @param m the data matrix defining this BiclusterExperiment. Should have rows
#'   as samples and features as columns
#' @return an instance of BiclusterExperiment-class
#'   containing the following slots, accessed with @@: 
#'   data: Object of class \code{\link{matrix}}. The original data. 
#'   annot: Object of class \code{\link{data.frame}}. Annotations provided by the user
#'   strategies: A \code{\link{list}} of \code{BiclusterStrategy} objects
#' @export
setGeneric("BiclusterExperiment", function(m, bcs = list(), phenoData = Biobase::annotatedDataFrameFrom(m, byrow = TRUE), featureData = annotatedDataFrameFrom(m, byrow = FALSE), pp = FALSE, maxNa = 0.5) {
  standardGeneric("BiclusterExperiment")
})

#' Careful, for m, rows are samples and columns are features
#' eSet objects store assayData transposed: rows are features and columns are samples.
#' For this reason I wrote a getter that returns a matrix with rows as samples, columns as features.
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
  ad$abund <- t(m)
  
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
#' Due to requirements of various biclustering methods, this function may warn
#' that the chosen biclustering method was overridden.
#' 
#' Also, if any matrix elements of abund(BiclusterExperiment) are missing, then
#' the BiclusterExperiment returned is not guaranteed to have the same number of
#' samples and features.
#' 
#' Optional boolean argument \code{duplicable} can force Spectral, Plaid, 
#' NIPALS, and SNMF to be duplicable. Alternatively, can allow ALS-NMF to be 
#' non-duplicable. Other arguments to specific functions in bicluster.R can be
#' provided.
#' 
#' @export
setGeneric("addStrat", function(bce, k, 
                                method,
                                maxNa = 1, duplicable = TRUE, silent = FALSE, 
                                ...) {
  standardGeneric("addStrat")
})
setMethod("addStrat", c(bce = "BiclusterExperiment", k = "numeric", 
                        method = "character"), 
          function(bce, k, method = c("als-nmf", "svd-pca", "snmf",
                                      "nipals-pca", "plaid", "spectral"), maxNa, 
                   duplicable, silent, ...) {
            # Validate parameters
            # k must be whole number, smaller than both dimensions of m
            m <- t(as.matrix(bce))
            k <- validateKM(k, m, method)
            if(!(maxNa <= 1 && maxNa >= 0)) {
              stop("Arg \"maxNa\" must be in the range of 0 to 1.")
            }
            
            if(method == "nipals-pca") {
              # to my knowledge, this is the only method so far that needs cleaning
              bce <- clean(bce, maxNa)
            }
            
            if(length(method) > 1) {
              stop("Argument \"method\" must be a single string")
            }
            method <- match.arg(method)
            # do this so that the recursive calls with a NULL method are also silent
            method.orig <- method
            
            if (method == "nipals-pca") {# Special code for NIPALS
              # Careful! because we're in tryCatch, warnings won't be printed.
              bcs <- tryCatch({
                BiclusterStrategy(m = m, k = k, method = method, ...)
              }, error = function(e) {
                if(method == "nipals-pca" &&
                   grepl(pattern = paste0("replacement has length zero"), 
                         x = e)) {
                  maxNa <- maxNa - (maxNa / 2)
                  message(paste("Cleaning with maxNAs at", maxNa))
                  
                  # Call recursively until success.
                  addStrat(bce, k, method.orig, maxNa, ...)
                } else { stop(e) }
              })
            } else {
              # Here, warnings and errors are thrown, not handled
              #### DEFUALT CALL ####
              bcs <- BiclusterStrategy(m = m, k = k, method = method, 
                                       duplicable = duplicable, ...)
            }
            
            name <- name(bcs)
            bce@strategies[[name]] <- bcs
            if(validObject(bce)) {
              message(paste("Added BiclusterStrategy named", name))
              return(bce)
            }
          })

# FIXME adapt from ExpressionSet method exprs so environment-style assayData can be accessed?
#' @export
setMethod("as.matrix", "BiclusterExperiment", function(x) {
  Biobase::assayDataElement(x, "abund")
})

setMethod("clean", c(object = "BiclusterExperiment"), function(object, maxNa) {
  results <- clean(t(as.matrix(object)), maxNa, TRUE)
  # [[2]] contains a vector of indexes of the remaining columns
  # [[1]] contains the cleaned matrix itself
  bce <- object[results[[2]][[2]], results[[2]][[1]]]
  strategies(bce) <- bce@strategies
  
  if(validObject(bce, test = FALSE)) return(bce)
})

#' Access a BiclusterStrategy contained by a BiclusterExperiment
#' 
#' Returns one BiclusterStrategy. Alternatively, a list of BiclusterStrategy
#' objects contained by a BiclusterExperiment can be retrieved using
#' \code{strategies(bce)}.
#' 
#' @param bce A BiclusterExperiment to access
#' @param id Either the integer index or the name of the BiclusterStrategy to 
#'   get
#' @export
setGeneric("getStrat", signature = "bce", function(bce, id) {standardGeneric("getStrat")})
setMethod("getStrat", c(bce = "BiclusterExperiment"), function(bce, id) {
  bce@strategies[[id]] 
})

#' Names of BiclusterStrategies in this BiclusterExperiment
#' @export
setMethod("names", "BiclusterExperiment", function(x) names(x@strategies))

#' Remove all BiclusterStrategy objects from this BiclusterExperiment
#' 
#' Returns a BiclusterExperiment containing only abundance data
#' 
#' @export
setGeneric("wipe", signature = "bce", function(bce) {standardGeneric("wipe")})
setMethod("wipe", c(bce = "BiclusterExperiment"), function(bce) {
  bce@strategies <- list()
  return(bce)
})

#' Removes all BiclusterStrategy objects except one
#' 
#' Returns a BiclusterExperiment encapsulating only the named BiclusterStrategy.
#' Can be useful to store the final results.
#' 
#' @export
setGeneric("wipeExcept", signature = c("bce", "bcs"),
           function(bce, bcs) {standardGeneric("wipeExcept")})
setMethod("wipeExcept", c(bce = "BiclusterExperiment", bcs = "numeric"),
           function(bce, bcs) {
             wipeExcept(bce, getStrat(bce, bcs))
           }
)
setMethod("wipeExcept", c(bce = "BiclusterExperiment", bcs = "character"),
           function(bce, bcs) {
             wipeExcept(bce, getStrat(bce, bcs))
           }
)
setMethod("wipeExcept", c(bce = "BiclusterExperiment",
                           bcs = "BiclusterStrategy"),
           function(bce, bcs) {
             if(any(names(bce) == name(bcs))) {
              bce@strategies <- list()
              bce@strategies[[name(bcs)]] <- bcs
              return(bce)
             } else {
               stop(paste("The given BiclusterStrategy is not contained by the",
"given BiclusterExperiment. Please create the BiclusterStrategy first."))
             }
           }
)
  
#' @export
setGeneric("strategies", signature = "bce", function(bce) {
  standardGeneric("strategies")
})
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
