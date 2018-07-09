setGeneric("as.matrix")

#### clean ####
#' Cleans a matrix by removing NAs
#'
#' @export
setGeneric("clean", function(object, maxNa = 0, dimsRemain = FALSE) {
  if(!(maxNa <= 1 && maxNa >= 0)) {
    stop("Arg \"maxNa\" must be in the range of 0 to 1.")
  }
  standardGeneric("clean")
})

#' Distance matrix accessor
#' 
#' Returns distance matrix as a class "matrix" object
#' 
#' @export
setGeneric("distance", signature = "bce", function(bce, ...) {standardGeneric("distance")})

#' @param bce A BiclusterExperiment to access
#' @param id Either the integer index or the name of the BiclusterStrategy to 
#'   get
#' @export
setGeneric("getStrat", signature = "bce", function(bce, id) {standardGeneric("getStrat")})

#' @export
setGeneric("name", signature = "bcs", function(bcs) {standardGeneric("name")})

#' Extract names
setGeneric("names")

#' Number of clusters in a BiclusterStrategy
#' 
#' @export
setGeneric("nclust", signature = "bcs", function(bcs) {standardGeneric("nclust")})

setGeneric("plot")

#' Bicluster prediction accessor
#' 
#' Get a binary matrix coding for bicluster membership. Rows represent samples and columns represent biclusters.
#' @export
setGeneric("pred", signature = "bcs", function(bcs) {standardGeneric("pred")})
