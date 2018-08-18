setGeneric("as.matrix")

#' Bicluster names
#' 
#' Get the names of all biclusters in a BiclusterStrategy. This function is
#' mostly useful for graphical functions.
setGeneric("bcNames", signature = "bcs", function(bcs, allBc = FALSE)
  standardGeneric("bcNames"))

setGeneric("biclust", function(bcs) standardGeneric("biclust"))
#' Clean a matrix or BiclusterExperiment
#'
#' Returns an object of the same type as \code{object}, containing all columns 
#' in \code{object} for which the fraction of present values is greater than
#' than \code{cleanParam}
#' 
#' @param object a matrix or BiclusterExperiment object
#' @param cleanParam a numeric within 0 to 1
#' @param dimsRemain if \code{TRUE}, returns a list of rows and columns
#'   remaining after cleaning. If \code{FALSE}, returns the cleaned object 
#'   itself.
#'
#' @export
setGeneric("clean", function(object, cleanParam = 0, dimsRemain = FALSE)
  standardGeneric("clean"))

#' Bicluster-feature matrix
#' 
#' Gets the logical matrix where rows represent features and columns represent
#' biclusters.
setGeneric("clusteredFeatures", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("clusteredFeatures")) 

#' Sample-bicluster matrix
#' 
#' Gets the logical matrix where rows represent biclusters and columns represent
#' samples.
setGeneric("clusteredSamples", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("clusteredSamples"))

#' Get a BiclusterStrategy
#' 
#' Returns one BiclusterStrategy in \code{bce} identified by \code{id}.
#' 
#' @param bce A \code{\linkS4class{BiclusterExperiment}} that contains one or
#'   more \code{\linkS4class{BiclusterStrategy}} objects
#' @param id A character string or numeric index identifying the desired
#'   \code{\linkS4class{BiclusterStrategy}}
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' addStrat(bce, k = 2)
#' getStrat(bce, "ALS-NMF | Otsu | 2")
setGeneric("getStrat", signature = "bce",
           function(bce, id) standardGeneric("getStrat"))

#' The biclustering algorithm used 
#'
#' Returns the name of the biclustering algorithm used to calculate a BiclusterStrategy
#' @export
setGeneric("method", signature = "bcs",
           function(bcs) standardGeneric("method"))

#' Loading matrix
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is k x N.
#'
#' @export
setGeneric("sampleFactor", signature = "bcs",
           function(bcs, allBc = FALSE)  standardGeneric("sampleFactor"))

#' Loading thresholds
#'
#' @export
setGeneric("loadingThresh", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("loadingThresh"))

#' Name a BiclusterStrategy
#' 
#' Gets a display-friendly name of a BiclusterStrategy that includes its
#' biclustering algorithm, its thresholding algorithm, and the number of
#' biclusters.
#' 
#' @param bcs a \code{\linkS4class{BiclusterStrategy}} object
setGeneric("name", signature = "bcs", function(bcs) standardGeneric("name"))

setGeneric("names")

#' @export
setGeneric("nclust", signature = "bcs", function(bcs, allBc = FALSE) 
  standardGeneric("nclust"))

setGeneric("plot")

#' Score matrix
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is M x k.
#''
#' @export
setGeneric("featureFactor", signature = "bcs", 
           function(bcs, allBc = FALSE) standardGeneric("featureFactor"))

#' Score thresholds
#'
#' @export
setGeneric("scoreThresh", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("scoreThresh"))

setGeneric("strategies", signature = "bce",
           function(bce) standardGeneric("strategies"))

setGeneric("strategies<-", 
           function(object, value) standardGeneric("strategies<-"))

setGeneric("threshAlgo", signature = "bcs",
           function(bcs) standardGeneric("threshAlgo"))

setGeneric("threshold", signature = c("m", "th"), 
           function(m, th, MARGIN = 2) standardGeneric("threshold"))

#' Wipe biclustering results
#'
#' Returns a BiclusterExperiment with all BiclusterStrategies removed
setGeneric("wipe", signature = "bce", function(bce) standardGeneric("wipe"))

setGeneric("wipeExcept", signature = c("bce"),
           function(bce, bcs) standardGeneric("wipeExcept"))
