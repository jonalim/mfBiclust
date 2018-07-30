setGeneric("as.matrix")

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
setGeneric("clean", function(object, cleanParam = 0, dimsRemain = FALSE) {
  standardGeneric("clean")
})

setGeneric("clusteredFeatures", signature = "bcs", function(bcs) {
  standardGeneric("clusteredFeatures")
}) 

setGeneric("clusteredSamples", signature = "bcs", function(bcs) {
  standardGeneric("clusteredSamples")
})

setGeneric("getStrat", signature = "bce", function(bce, id) {standardGeneric("getStrat")})

#' The biclustering algorithm used 
#'
#' Returns the name of the biclustering algorithm used to calculate a BiclusterStrategy
#' @export
setGeneric("method", signature = "bcs", function(bcs) {
  standardGeneric("method")
})

#' Loading matrix
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is k x N.
#'
#' @export
setGeneric("loading", signature = "bcs", function(bcs) {
  standardGeneric("loading")
})

#' Loading thresholds
#'
#' @export
setGeneric("loadingThresh", signature = "bcs", function(bcs) {
  standardGeneric("loadingThresh")
})

setGeneric("name", signature = "bcs", function(bcs) {standardGeneric("name")})

setGeneric("names")

#' @export
setGeneric("nclust", signature = "bcs", function(bcs) {standardGeneric("nclust")})

setGeneric("plot")

#' Score matrix
#'
#' For a data matrix M x N factorized to produce k biclusters, the score matrix is M x k.
#''
#' @export
setGeneric("score", signature = "bcs", function(bcs) {
  standardGeneric("score")
})

#' Score thresholds
#'
#' @export
setGeneric("scoreThresh", signature = "bcs", function(bcs) {
  standardGeneric("scoreThresh")
})

setGeneric("strategies", signature = "bce", function(bce) {
  standardGeneric("strategies")
})

setGeneric("strategies<-", function(object, value) standardGeneric("strategies<-"))

setGeneric("threshAlgo", signature = "bcs", function(bcs) {
  standardGeneric("threshAlgo")
})

setGeneric("threshold", signature = c("m", "th"), function(m, th, MARGIN = 2) {
  standardGeneric("threshold")
})

setGeneric("wipe", signature = "bce", function(bce) {standardGeneric("wipe")})

setGeneric("wipeExcept", signature = c("bce"),
           function(bce, bcs) {standardGeneric("wipeExcept")})
