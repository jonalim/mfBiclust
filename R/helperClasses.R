#' A pair of matrix factors
#' 
#' Two matrices that can be thresholded to determine which samples and
#' features are members of which biclusters. After running a
#' matrix-factorization biclustering method, such as SVD-PCA, ALS-NMF, or SNMF,
#' these matrices will be score and loading matrices.
#' @slot W an m x k matrix
#' @slot H a k x n matrix
#' @export
setClass("genericFactorization", slots = list(W = "matrix", H = "matrix"))

#' A biclustering algorithm result
#' 
#' Contains the results and metadata from running a biclustering algorithm.
#' 
#' @slot fit a \code{\link{genericFactorization-class}} object
#' @slot method the name of the biclustering algorithm used
#' @export
setClass("genericFit", slots = list(fit = "genericFactorization", 
                                    method = "character"))