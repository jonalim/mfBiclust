#### clean ####
#' Clean a matrix
#'
#' Returns an object of the same type as \code{object}, containing all columns 
#' in \code{object} for which the fraction of missing values is less than
#' \code{cleanParam}
#'
#' @export
setGeneric("clean", function(object, cleanParam = 0, dimsRemain = FALSE) {
  standardGeneric("clean")
})

#' @export
setGeneric("nclust", signature = "bcs", function(bcs) {standardGeneric("nclust")})

setGeneric("plot")

