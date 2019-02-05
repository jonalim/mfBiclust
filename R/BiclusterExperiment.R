#' @include BiclusterStrategy.R
#' @include generics.R
NULL

#### CLASS #####################################################################
#' Class "BiclusterExperiment" for data and biclusters
#'
#' This class encapsulates data for one or more biclustering runs derived from
#' the same abundance data. Objects can be created using the
#' \code{\link{BiclusterExperiment}} constructor. A subclass of
#' \code{\link{eSet}}.
#'
#' @slot assayData Object of class \code{\link{AssayData-class}} in the
#'   "list" storage mode. \code{AssayData} must have a matrix named
#'   "abund" with rows representing features and columns representing samples.
#'   Any other matrices in assayData are ignored by this package.
#' @slot strategies A \code{\link{list}} of
#'   \code{\link{BiclusterStrategy-class}} objects
#' @slot phenoData See \code{\link{eSet-class}}.
#' @slot featureData \code{\link{eSet-class}}.
#' @slot experimentData \code{\link{eSet-class}}. Not accessed by any
#'   functions in \code{\link{mfBiclust}}.
#' @slot annotation \code{\link{eSet-class}}. Not accessed by any
#'   functions in \code{\link{mfBiclust}}.
#' @slot protocolData \code{\link{eSet-class}}. Not accessed by any
#'   functions in \code{\link{mfBiclust}}.
#'
#' @param x A \code{BiclusterExperiment-class} object
#'
#' @seealso \code{\link{eSet-class}}
#' @seealso \code{\link{getStrat}}
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#'
#' @importClassesFrom Biobase eSet
setClass("BiclusterExperiment", slots = list(
  strategies = "list"
), contains = "eSet")

# Coercing ExpressionSet (gene expression) to BiclusterExperiment
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

# Only possible after running groupChromPeaks...which should be run after
# adjusting retention times.
setAs("XCMSnExp", "BiclusterExperiment", function(from) {
  if(!requireNamespace("xcms", quietly = TRUE)) {
    stop("Please install XCMS")
  }
  ad <- xcms::featureValues(from)
  pd <- xcms::phenoData(from)
  # phenodata sample names are set to column names of featureValues(from)
  if(any(rownames(pd@data) != colnames(ad))) {
    rownames(pd@data) <- colnames(ad)
  }
  bce <- BiclusterExperiment(m = ad,
             phenoData = pd,
             featureData = as.data.frame(xcms::featureDefinitions(from)))
  bce@experimentData <- experimentData(from)
  # missing from@protocolData, but I don't know how to make that object
  # compatible with other objects in an eSet
  if(validObject(bce, test = FALSE)) bce
})

#### CONSTRUCTOR ###############################################################
#' Construct a container for biclustering runs
#'
#' Constructs a \code{\link{BiclusterExperiment-class}} object holding data,
#' with a slot for future biclustering runs performed on that data. To run
#' biclustering, use \code{\link{addStrat}()}.
#'
#' @param m The data on which derived BiclusterStrategy objects will be
#'   calculated. Should be a matrix or a class coercible to matrix, with samples
#'   as columns and features as rows
#' @param bcs May be a \code{list} of BiclusterStrategy objects that have been
#'   computed from \code{m} previously
#' @param phenoData Metadata about the rows of \code{m}
#' @param featureData Metadata about the columns of \code{m}
#' @param pp Whether to preprocess \code{m} by removing rows and columns with
#'   excessive missing data
#' @param maxNa If \code{pp}, rows and columns with over \code{maxNa} elements
#'   missing are discarded
#'
#' @return A \code{\link{BiclusterExperiment-class}} object containing the data
#'  passed to the constructor
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#'
#' @export
setGeneric(
    "BiclusterExperiment",
    function(m,
             bcs = list(),
             phenoData = Biobase::annotatedDataFrameFrom(m, byrow = FALSE),
             featureData = annotatedDataFrameFrom(m, byrow = TRUE),
             pp = FALSE,
             maxNa = 0) {
    standardGeneric("BiclusterExperiment")
})

#' @inheritParams BiclusterExperiment
#' @describeIn BiclusterExperiment Default constructor method
setMethod(
    "BiclusterExperiment", c(m = "matrix"),
    function(m, bcs, phenoData, featureData, pp, maxNa) {
        if(pp) {
            m <- clean(m, 1 - maxNa)
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

#' @inheritParams BiclusterExperiment
#' @describeIn BiclusterExperiment Attempts to coerce \code{m} to a matrix before
#'  encapsulating it in a \code{BiclusterExperiment}
setMethod("BiclusterExperiment", c(m = "ANY"), function(m, bcs, phenoData, featureData, pp, maxNa) {
  BiclusterExperiment(as.matrix(m), bcs, phenoData, featureData, pp, maxNa)
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
    if(!all(names(object@strategies) == vapply(object@strategies, name, character(1)))) {
      msg <- c(msg, paste("List names of object@strategies are not identical",
                          "to BiclusterStrategy names"))
    } else {
      validBcs <- vapply(names(object), FUN.VALUE = logical(1), FUN = function(bcs) {
        inherits(getStrat(object, bcs), "BiclusterStrategy")
      })
      if (!all(validBcs)) {
        msg <- c(msg, "All strategies must be BiclusterStrategy objects.")
      } else {
        newMsgs <- lapply(names(object), FUN = function(bcs) { # Check validity of all strategies
          res <- validObject(getStrat(object, bcs), test = TRUE)
          if(inherits(res, "character")) { return(paste0(res, "\n")) } else { return(NULL) }
        })
        msg <- c(msg, unlist(newMsgs))
      }
    }
  } else {
    msg <- c(msg, "The strategies slot must be a 'list' object")
  }
  if(is.null(msg)) TRUE else msg
}
setValidity("BiclusterExperiment", validBiclusterExperiment)

#' @rdname BiclusterExperiment-class
#' @export
setMethod("as.matrix", "BiclusterExperiment", function(x) {
  Biobase::assayDataElement(x, "abund")
})

#' @rdname clean
setMethod("clean", c(object = "BiclusterExperiment"),
          function(object, cleanParam, dimsRemain = FALSE) {
  if(!(cleanParam <= 1 && cleanParam >= 0)) {
    stop("Arg \"cleanParam\" must be in the range of 0 to 1.")
  }
  results <- clean(as.matrix(object), cleanParam, TRUE)
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

#' @rdname getStrat
#' @export
setMethod("getStrat", c(bce = "BiclusterExperiment"), function(bce, id) {
  res <- strategies(bce)[[id]]
  if(is.null(res)) { #Occurs if is.character(id) and no matching BCS is found
    stop("BiclusterStrategy not found")
  }
  return(res)
})

#' @inheritParams BiclusterExperiment-class
#' @describeIn BiclusterExperiment Retrieves the names of any encapsulated
#'  \code{\link{BiclusterStrategy-class}} objects.
#' @export
setMethod("names", "BiclusterExperiment", function(x) names(x@strategies))

#' @inheritParams BiclusterExperiment-class
#' @describeIn BiclusterExperiment Get/set a list of the BiclusterStrategy objects
#'   contained by this BiclusterExperiment.
#' @export
setMethod("strategies", c(x = "BiclusterExperiment"), function(x) {
  x@strategies
})

setReplaceMethod("strategies", signature(object = "BiclusterExperiment",
                                         value = "list"),
                 function(object, value) {
                   object@strategies <- value
                   if(validObject(object, test = FALSE)) {
                     object }
                 }
)

#' @rdname wipe
#' @export
setMethod("wipe", c(bce = "BiclusterExperiment"), function(bce) {
  strategies(bce) <- list()
  return(bce)
})

#' @rdname wipeExcept
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
