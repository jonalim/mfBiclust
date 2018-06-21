setGeneric("as.matrix")

#### clean ####
#' Cleans a matrix or BiclusterExperiment
#'
#' @export
setGeneric("clean", function(object, maxNa = 0, dimsRemain = FALSE) {
  if(!(maxNa <= 1 && maxNa >= 0)) {
    stop("Arg \"maxNa\" must be in the range of 0 to 1.")
  }
  standardGeneric("clean")
})
setMethod("clean", c(object = "matrix"), function(object, maxNa, dimsRemain) {
  maxNaPerRow <- round(maxNa * ncol(object))
  maxNaPerCol <- round(maxNa * nrow(object))
  
  # Both 0s and NAs can foul up NIPALS
  goodRows <- apply(object, MARGIN = 1, function(row) 
    (sum(is.na(row)) + sum(row == 0, na.rm = TRUE)) < maxNaPerRow)
  goodCols <- apply(object, MARGIN = 2, function(col) 
    (sum(is.na(col)) + sum(col == 0, na.rm = TRUE)) < maxNaPerCol)
  
  object <- object[goodRows, goodCols]
  
  if(dimsRemain) {
    list(obj = object, dimsRemain = list(goodRows, goodCols))
  } else { object }
  
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
