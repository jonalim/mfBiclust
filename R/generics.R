setGeneric("as.matrix")

#' Names of biclusters in this BiclusterStrategy
#'
#' Returns a character vector of bicluster names, which can be used to reference
#' columns of \code{featureFactor(bcs)} and \code{clusteredFeatures(bcs)}, and
#' rows of \code{sampleFactor(bcs)} and code{clusteredSamples(bcs)}. Mostly
#' useful for graphical functions.
#'
#' @param bcs A BiclusterStrategy
#' @param allBc TRUE to return data from biclusters that were not requested
#'   during the original call to \code{\link{addStrat}}
#'
#' @return A character vector.
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' bcNames(bcs)
setGeneric("bcNames", signature = "bcs", function(bcs, allBc = FALSE)
  standardGeneric("bcNames"))

setGeneric("biclust", function(bcs) standardGeneric("biclust"))

#' Clean a matrix or BiclusterExperiment
#'
#' Returns an object of the same type as \code{object}, containing all columns
#' in \code{object} for which the fraction of present values is greater than
#' than \code{cleanParam}
#'
#' @param object a numeric \code{matrix} or
#'   \code{\link{BiclusterExperiment-class}} object
#' @param cleanParam a numeric within 0 to 1
#' @param dimsRemain if \code{TRUE}, returns a list of rows and columns
#'   remaining after cleaning. If \code{FALSE}, returns the cleaned object
#'   itself.
#'
#' @return A \code{BiclusterExperiment-class} or a \code{matrix}
setGeneric("clean", function(object, cleanParam = 0, dimsRemain = FALSE)
  standardGeneric("clean"))

#' Bicluster-feature matrix
#'
#' Returns a logical matrix \eqn{A_{m x k}} where element \eqn{A_{i,j}} is TRUE
#' if feature \eqn{i} is in bicluster \eqn{j}. Elsewhere in \code{mfBiclust}
#' documentation, this format is called a "membership matrix".
#'
#' @param bcs A BiclusterStrategy
#' @param allBc TRUE to return data from biclusters that were not requested
#'   during the original call to \code{\link{addStrat}}
#'
#' @return A logical matrix
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' clusteredFeatures(bcs)
#'
setGeneric("clusteredFeatures", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("clusteredFeatures"))

#' Sample-bicluster matrix
#'
#' Returns a logical matrix \eqn{B_{k x n}} where element \eqn{B_{i,j}} is TRUE
#' if sample \eqn{j} is in bicluster \eqn{i}. Elsewhere in \code{mfBiclust}
#' documentation, this format is called a "membership matrix".
#'
#' @param bcs A BiclusterStrategy
#' @param allBc TRUE to return data from biclusters that were not requested
#'   during the original call to \code{\link{addStrat}}
#'
#' @return A logical matrix
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' clusteredSamples(bcs)
#'
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
#'
#' @return A \code{\link{BiclusterStrategy-class}} object
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2)
#' getStrat(bce, "ALS-NMF | Otsu | 2")
#' getStrat(bce, 1)
setGeneric("getStrat", signature = "bce",
           function(bce, id) standardGeneric("getStrat"))

#' The biclustering algorithm used
#'
#' Returns the name of the biclustering algorithm used to calculate a
#' BiclusterStrategy
#' @export
#'
#' @param bcs A \code{\link{BiclusterStrategy-class}} object
#'
#' @return character string
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2)
#' bcs <- getStrat(bce, 1)
#' method(bcs)
setGeneric("method", signature = "bcs",
           function(bcs) standardGeneric("method"))

#' Sample matrix factor
#'
#' For a data matrix \eqn{A_{m x n}} factorized as \eqn{W_{m x k}H_{k x n}},
#' returns \eqn{H}.
#'
#' @param bcs A Bic
#'
#' @return a numeric matrix.
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' sampleFactor(bcs)
#'
#' @export
setGeneric("sampleFactor", signature = "bcs",
           function(bcs, allBc = FALSE)  standardGeneric("sampleFactor"))

#' Sample thresholds
#'
#' Returns the thresholds that were applied to each row of the factor matrix
#' \eqn{H_{k x n}} to determine the samples in each bicluster.
#'
#' @return a numeric vector of length \eqn{k}
#'
#'
#' @param bcs A BiclusterStrategy
#' @param allBc TRUE to return data from biclusters that were not requested
#'   during the original call to \code{\link{addStrat}}
#'
#' @export
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' sampleThresh(bcs)
setGeneric("sampleThresh", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("sampleThresh"))


#' Name of a BiclusterStrategy
#'
#' Gets a display-friendly name of a BiclusterStrategy that includes its
#' biclustering algorithm, its thresholding algorithm, and the number of
#' biclusters.
#'
#' This function may not be used to modify a BiclusterStrategy's name.
#' @export
#'
#' @param bcs A \code{\link{BiclusterStrategy-class}} object
#'
#' @return character string
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2)
#' bcs <- getStrat(bce, 1)
#' name(bcs)
#' @rdname name
setGeneric("name", signature = "bcs", function(bcs) standardGeneric("name"))

setGeneric("names")

#' Number of biclusters
#'
#' Returns the number of biclusters in a \code{\link{BiclusterStrategy-class}}
#' object
#'
#' @param bcs A BiclusterStrategy
#' @param allBc TRUE to return the total number of biclusters that were found
#'   during the original call to \code{\link{addStrat}}
#'
#' @return integer
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2)
#' bcs <- getStrat(bce, 1)
#' nclust(bcs)
#'
#' @export
setGeneric("nclust", signature = "bcs", function(bcs, allBc = FALSE)
  standardGeneric("nclust"))

#' Expression plot
#'
#' Create a heatmap from assay data in a BiclusterExperiment object. A BiclusterStrategy and bicluster must be named if ordering = "cluster" or
#' annotateCluster = TRUE.
#'
#' @param x A \code{link{BiclusterExperiment-class}} object
#' @param y Not used in the \code{BiclusterStrategy} method
#' @param ... Additional arguments:\describe{
#'   \item{logBase}{If TRUE, every assay value is transformed by the natural log}
#'   \item{phenoLabels}{Display phenotype labels along the samples?}
#'   \item{biclustLabels}{Display bicluster labels along the samples and features?}
#'   \item{ordering}{If "input", samples and features are not reordered. If
#'   "distance", sampels and features are reordered based on Euclidean distance.
#'   If "cluster", samples and features are reordered based on bicluster
#'   membership.}
#'   \item{strategy}{The BiclusterStrategy to read biclusters from, if bicluster-
#'  based information is to be plotted.}
#'   \item{rowNames}{TRUE to show feature names}
#'   \item{colNames}{TRUE to show sample names}
#' }
#'
#' @return a \code{\link[pheatmap]{pheatmap}-class} object
#'
#' @examples
#' bce <- BiclusterExperiment(cancer_benchmark[[1]]$data,
#'                            phenoData = cancer_benchmark[[1]]$labels)
#' phenoData(bce)
#' plot(bce, phenoLabels = "phenotype")
setGeneric("plot")

#' Feature matrix factor
#'
#' For a data matrix \eqn{A_{m x n}} factorized as \eqn{W_{m x k}H_{k x n}},
#' returns \eqn{W}.
#'
#' @param bcs A BiclusterStrategy
#' @param allBc TRUE to return data from biclusters that were not requested
#'   during the original call to \code{\link{addStrat}}
#'
#' @return a numeric matrix
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' featureFactor(bcs)
#'
#' @export
setGeneric("featureFactor", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("featureFactor"))

#' Feature thresholds
#'
#' Returns the thresholds that were applied to each row of the feature matrix
#' \eqn{W_{m x k}} to determine the samples in each bicluster.
#'
#' @param bcs A BiclusterStrategy
#' @param allBc TRUE to return data from biclusters that were not requested
#'   during the original call to \code{\link{addStrat}}
#'
#' @return a numeric vector of length \eqn{m}
#'
#' @export
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' featureThresh(bcs)
setGeneric("featureThresh", signature = "bcs",
           function(bcs, allBc = FALSE) standardGeneric("featureThresh"))

#' Results of biclustering runs on a single dataset
#'
#' Returns the contents of \code{\link{BiclusterExperiment-class}@@strategies}.
#'
#' @param bce A \code{BiclusterExperiment-class} object
#'
#' @value A list containing zero or more \code{\link{BiclusterStrategy-class}}
#'   objects
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bce <- addStrat(bce, k = 3, method = "als-nmf")
#' strategies(bce)
setGeneric("strategies", signature = "bce",
           function(bce) standardGeneric("strategies"))

setGeneric("strategies<-",
           function(object, value) standardGeneric("strategies<-"))

#' Algorithm used to calculate a threshold
#'
#' Returns \code{\link{BiclusterStrategy-class}@@threshAlgo}
#'
#' @param bcs A \code{BiclusterStrategy-class} object
#'
#' @return A character string
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' threshAlgo(bcs)
setGeneric("threshAlgo", signature = "bcs",
           function(bcs) standardGeneric("threshAlgo"))

#' Apply thresholds to a score or loading matrix
#'
#' Returns a binary matrix of the same size as \code{m} where all elements over
#' the threshold are 1.
#'
#' @param m A numeric matrix
#' @param th A numeric vector with the same number of elements as
#'   \code{dim(m)[[MARGIN]]}
#' @param MARGIN 1 to threshold rows; 2 to threshold columns.
#'
#' @return A binary numeric matrix of the same dimensions as \code{m}
#'
#' @examples
#' m <- matrix(rnorm(100), 10, 10)
#' m
#' threshold(m = m, th = rep(0, 10), MARGIN = 1)
#'
#' @export
setGeneric("threshold", signature = c("m", "th"),
           function(m, th, MARGIN = 2) standardGeneric("threshold"))

#' Wipe biclustering results
#'
#' Returns a \code{\link{BiclusterExperiment-class}} object with all
#' BiclusterStrategies removed.
#'
#' @param bce A \code{BiclusterExperiment-class} object
#'
#' @return A \code{BiclusterExperiment-class} object
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' wipe(bce)
setGeneric("wipe", signature = "bce", function(bce) standardGeneric("wipe"))

#' Wipe all biclustering results except one
#'
#' Returns a BiclusterExperiment with the same data and metadata as \code{bce},
#' but with \code{BiclusterExperiment@strategies} containing the single
#' specified \code{BiclusterStrategy} instance. The \code{BiclusterStrategy} to
#' preserve can be specified by name, by index, or by reference to an identical
#' \code{BiclusterStrategy}. This function is useful for isolating your
#' "favorite" biclustering result.
#'
#' @param bce A \code{BiclusterExperiment-class} object
#' @param bcs A \code{BiclusterStrategy-class} object in \code{bce@@strategies},
#'   or the name or index of an object in \code{bce@@strategies}
#'
#' @return A BiclusterExperiment object
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bce <- addStrat(bce, k = 5, method = "als-nmf")
#' wipeExcept(bce, 2) # preserves the second BiclusterStrategy
setGeneric("wipeExcept", signature = c("bce"),
           function(bce, bcs) standardGeneric("wipeExcept"))
