#' @include biclusterGUI.R
#' @include BiclusterExperiment.R
#' @include helperFunctions.R
NULL

#' Analyze a matrix and display results graphically
#'
#' Runs all steps of the biclustering pipeline.
#'
#' \code{mfbc} is the fastest way to analyze data. The only required option is
#' the input matrix; the default pipeline uses non-negative matrix factorization
#' with otsu thresholding. To compare multiple pipelines, manual construction of
#' a \code{\link{BiclusterExperiment}} object is recommended.
#'
#' @param x The input matrix
#' @param annotations A dataframe, where each column contains a set of labels on
#'   the rows of \code{x}
#' @param ks the number of clusters to create, can be atomic or a vector
#' @param bicluster the biclustering algorithm
#' @param scoreThresh the thresholding algorithm to use for scores. Ignored if
#'   \code{bicluster} is \code{"plaid"} or \code{"bimax"}. Currently "snmf/l"
#'   and "pca" are implemented
#' @param loadingThresh the thresholding algorithm to use for loading. See note
#'   on \code{scoreThresh}
#' @param bcv Perform bi-cross-validation? Not yet implemented.
#'
#' @seealso \code{\link{BiclusterExperiment}}
#'
#' @export
setGeneric("mfbc", signature = "m", function(m, ...) {
  standardGeneric("mfbc")
})

setMethod(
  "mfbc",
  c(m = "matrix"),
  function(m, annotations = NULL, ks = 2:5, 
           bicluster = c("snmf/l", "pca"), 
           scoreThresh = c("otsu"),
           loadingThresh = c("otsu"), 
           bcv = FALSE, nogui = FALSE) {
    # Save old options; force warnings to appear interactively. (Does Shiny
    # require this to display interactive dialogs???)
    old <- options(stringsAsFactors = FALSE)
    on.exit(options(old), add = TRUE)
    options("warn" = 1)

    # invisible dummy matrix, if input is in wrong format.
    # TODO import cleaning for matrices
    # if (any(m < 0)) {
    #   m <- matrix(sample(x = 1:100, size = 600, replace = TRUE), nrow = 20)
    # }
    dataset <- as.matrix(m)

    # dummy annotations
    annotations <- data.frame(species = sample(x = c("human", "mouse"), 
                                               size = nrow(dataset), 
                                               replace = TRUE), 
                              tissue = sample(x = c("brain", "kidney", 
                                                    "cancer"), 
                                              size = nrow(dataset), 
                                              replace = TRUE),
                              tissue2 = sample(x = c("brain", "kidney", 
                                                    "cancer"), 
                                              size = nrow(dataset), 
                                              replace = TRUE),
                              tissue3 = sample(x = c("brain", "kidney", 
                                                    "cancer"), 
                                              size = nrow(dataset), 
                                              replace = TRUE),
                              species2 = sample(x = c("human", "mouse"), 
                                               size = nrow(dataset), 
                                               replace = TRUE),
                              species3 = sample(x = c("human", "mouse"), 
                                               size = nrow(dataset), 
                                               replace = TRUE))
  
    ks <- validateKs(ks)
    bicluster <- match.arg(bicluster)

    strats <- lapply(ks, function(k) BiclusterStrategy(m, k, bicluster, 
                                                              scoreThresh, 
                                                              loadingThresh))
    names(strats) <- as.character(ks)

    bce <- BiclusterExperiment(m = dataset, bcs = strats, annot = annotations, 
                               bcv = FALSE)
    if (bcv == TRUE) {
      warning("Bi-cross-validation is still under development. msNMF
        cannot predict the optimal clustering strategy.")
    }

    if(nogui) {
      bce
    } else {
      biclusterGUI(bce)
    }
  }
)
