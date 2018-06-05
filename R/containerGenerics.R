setGeneric("as.matrix")
setGeneric("plot")

#' Extract names
setGeneric("names")

#' @export
setGeneric("getStrat", signature = "bce", function(bce, ...) {standardGeneric("getStrat")})

#' @export
setGeneric("name", signature = "bcs", function(bcs) {standardGeneric("name")})

#' @export
setGeneric("nclust", signature = "bcs", function(bcs) {standardGeneric("nclust")})
#' nclust convenience method
#' 
#' Helper in case nclust is called on a list containing one BiclusterStrategy
setMethod("nclust", c(bcs = "list"), function(bcs) {
  ncol(bcs[[1]]@factors@fit@W)
}
)

#' Bicluster prediction accessor
#' 
#' Get a binary matrix coding for bicluster membership. Rows represent samples and columns represent biclusters.
#' @export
setGeneric("pred", signature = "bcs", function(bcs) {standardGeneric("pred")})

#' Distance matrix accessor
#' 
#' Returns distance matrix as a class "matrix" object
#' 
#' @export
setGeneric("distMat", signature = "bce", function(bce, ...) {standardGeneric("distMat")})
