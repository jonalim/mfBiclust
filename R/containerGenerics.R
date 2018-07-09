setGeneric("as.matrix")

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
