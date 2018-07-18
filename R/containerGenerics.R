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

#' @export
setGeneric("name", signature = "bcs", function(bcs) {standardGeneric("name")})

#' Extract names
setGeneric("names")

#' Number of clusters in a BiclusterStrategy
#' 
#' @export
setGeneric("nclust", signature = "bcs", function(bcs) {standardGeneric("nclust")})

setGeneric("plot")

