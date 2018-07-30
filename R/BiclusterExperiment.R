#' @include BiclusterStrategy.R
#' @include generics.R
NULL

#### CLASS #####################################################################
#' Class "BiclusterExperiment" for data and biclusters
#'
#' This class encapsulates data for one or more biclustering runs derived from
#' the same abundance data. Objects can be created using the
#' \code{\link{BiclusterExperiment}} constructor. A subclass of
#' \code{\link[Biobase]{eSet}}.
#'
#' @slot assayData Object of class \code{\link[Biobase]{AssayData-class}} in the
#'   "list" storage mode. \code{assayData} must have a matrix named
#'   "abund" with rows representing features and columns representing samples.
#'   Any other matrices in assayData are ignored by this package.
#' @slot strategies A \code{\link{list}} of
#'   \code{\link{BiclusterStrategy-class}} objects
#' @slot phenoData See \code{\link[Biobase]{eSet}}.
#' @slot featureData \code{\link[Biobase]{eSet}}.
#' @slot experimentData \code{\link[Biobase]{eSet}}. Not accessed by any 
#'   functions in \code{\link{mfBiclust}}.
#' @slot annotation \code{\link[Biobase]{eSet}}. Not accessed by any 
#'   functions in \code{\link{mfBiclust}}.
#' @slot protocolData \code{\link[Biobase]{eSet}}. Not accessed by any 
#'   functions in \code{\link{mfBiclust}}.
#'
#' @seealso \code{\link[Biobase]{eSet}}
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
  
  if(validObject(bce, test = FALSE)) {
    if(dimsRemain) {
      return(list(obj = bce, dimsRemain = results[[2]]))
    } else { bce }
    }
})

#' @describeIn BiclusterExperiment Get a BiclusterStrategy contained by a
#'   BiclusterExperiment, by providing either name or integer index
#' @export
setMethod("getStrat", c(bce = "BiclusterExperiment"), function(bce, id) {
  strategies(bce)[[id]] 
})

#' @describeIn BiclusterExperiment Character names of BiclusterStrategies in this BiclusterExperiment
#' @export
setMethod("names", "BiclusterExperiment", function(x) names(x@strategies))

#' @describeIn BiclusterExperiment Return this BiclusterExperiment with all BiclusterStrategy objects removed
#' @export
setMethod("wipe", c(bce = "BiclusterExperiment"), function(bce) {
  strategies(bce) <- list()
  return(bce)
})

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

#' @describeIn BiclusterExperiment Get/set a list of the BiclusterStrategy objects
#'   contained by this BiclusterExperiment.
#' @export
setMethod("strategies", c(bce = "BiclusterExperiment"), function(bce) {
  bce@strategies
})

setReplaceMethod("strategies", signature(object = "BiclusterExperiment",
                                         value = "list"),
                 function(object, value) { 
                   object@strategies <- value 
                   if(validObject(object, test = FALSE)) {
                     object }
                 }
)
