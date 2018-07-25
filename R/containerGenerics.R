setGeneric("as.matrix")

#### clean ####
#' Cleans a matrix by removing NAs
#'
#' @export
setGeneric("clean", function(object, cleanParam = 0, dimsRemain = FALSE) {
  standardGeneric("clean")
})

#' Distance matrix accessor
#' 
#' Returns distance matrix as a class "matrix" object
#' 
#' @export
setGeneric("distance", signature = "bce", function(bce, ...) {standardGeneric("distance")})

#' @export
setGeneric("name", signature = "bcs", function(bcs) {standardGeneric("name")})

#' Extract names
setGeneric("names")

#' Number of clusters in a BiclusterStrategy
#' 
#' @export
setGeneric("nclust", signature = "bcs", function(bcs) {standardGeneric("nclust")})

setGeneric("plot")

