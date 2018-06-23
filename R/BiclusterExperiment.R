#' @include BiclusterStrategy.R
#' @include helperFunctions.R
#' @include containerGenerics.R
NULL

#### CLASS #####################################################################
#' Class "BiclusterExperiment" for multiple biclustering results
#'
#' This class encapsulates factorization and thresholding data for one or more
#' biclustering runs. Objects can be created using the
#' \code{\link{BiclusterExperiment}} constructor.
#'
#' data = "matrix", annot = "data.frame", strategies = "list", distance = "matrix"
#'
#' @slot data Object of class \code{\link{matrix}}. The original data.
#' @slot annot Object of class \code{\link{data.frame}}. Annotations provided by the user
#' @slot strategies A \code{\link{list}} of \code{BiclusterStrategy} objects
#' @slot distance Object of class \code{\link{matrix}}. Pairwise distances
#'   between the rows of \code{data}
#' @importClassesFrom Biobase eSet
setClass("BiclusterExperiment", slots = list(
  strategies = "list", distance = "dist"
  ### FIXME: Change to similarity matrix? Use cor(t(as.matrix(x)), method = "pearson") to get similarity. Use "clustering_distance_rows = "correlation" in pheatmap calls.
), contains = "eSet")

setAs("ExpressionSet", "BiclusterExperiment", function(from) {
  ad <- Biobase::exprs(from)
  
  # remove genes with any NA
  # naIndex <- which(rowSums(is.na(ad)) > 0)
  # from <- from[-naIndex, ]
  # ad <- Biobase::exprs(from)
  
  d <- dist(t(Biobase::exprs(from)), method = "euclidean")
  
  # Add "abund" matrix and remove "exprs" matrix from the assayData object
  from <- Biobase::assayDataElementReplace(from, "abund", ad, validate = FALSE)
  from <- Biobase::assayDataElementReplace(from, "exprs", NULL)
  bce <- new("BiclusterExperiment", assayData = Biobase::assayData(from), 
             phenoData = Biobase::phenoData(from), 
             featureData = Biobase::featureData(from), strategies = list(), 
             distance = d)
  if(validObject(bce, test = FALSE)) bce
})

#### CONSTRUCTOR ###############################################################
#' Perform multiple biclustering runs
#'
#' BiclusterExperiment constructs an object holding data from multiple
#' biclustering runs.
#'
#' This function is useful for constructing one BiclusterExperiment
#' encapsulating results of different pipelines (e.g. comparing NMF with PCA).
#' For comparing results from the same pipeline with differing values of
#' \code{k}, \code{\link{mfbc()}} is recommended.
#' 
#' @param m the data matrix defining this BiclusterExperiment. Should have rows
#'   as samples and features as columns
#' @return an instance of BiclusterExperiment-class
#'   containing the following slots, accessed with @@: 
#'   data: Object of class \code{\link{matrix}}. The original data. 
#'   annot: Object of class \code{\link{data.frame}}. Annotations provided by the user
#'   strategies: A \code{\link{list}} of \code{BiclusterStrategy} objects
#'   distance: Object of class \code{\link{matrix}}. Pairwise distances
#'     between the rows of \code{data}
setGeneric("BiclusterExperiment", function(m, ...) {
  standardGeneric("BiclusterExperiment")
})

#' Careful, for m, rows are samples and columns are features
#' eSet objects store assayData transposed: rows are features and columns are samples.
#' For this reason I wrote a getter that returns a matrix with rows as samples, columns as features.
setMethod("BiclusterExperiment", c(m = "matrix"), function(m, bcs = list(), phenoData = Biobase::annotatedDataFrameFrom(m, byrow = TRUE), featureData = annotatedDataFrameFrom(m, byrow = FALSE), bcv = FALSE, maxNa = 0.5) {
  if (bcv == TRUE) {
    warning("Bi-cross-validation is still under development. msNMF
                      cannot predict the optimal clustering strategy.")
  }
  m <- clean(m, maxNa)
  d <- dist(m, method = "euclidean")
  
  if(!inherits(phenoData, "AnnotatedDataFrame")) {
    phenoData <- AnnotatedDataFrame(phenoData)
  }
  if(!inherits(featureData, "AnnotatedDataFrame")) {
    featureData <- AnnotatedDataFrame(featureData)
  }
  
  ad <- Biobase::assayDataNew(storage.mode = "list")
  ad$abund <- t(m)
  
  if(inherits(bcs, "list")) {
    names(bcs) <- lapply(bcs, function(bcs) {name(bcs)})
  } else if(!is.null(bcs)) {
    bcs <- list(bcs)
    names(bcs) <- name(bcs[[1]])
  } else {
    bcs <- list()
  }
  new("BiclusterExperiment", assayData = ad, phenoData = phenoData, featureData = featureData, strategies = bcs, distance = d)
})

#### METHODS ###################################################################

validBiclusterExperiment <- function( object ) {
  msg <- NULL
  if(!inherits(object, "BiclusterExperiment")) {
    msg <- c(msg, paste("Cannot validate a", class(object), 
                        "as BiclusterExperiment"))
  }
  if(!"abund" %in% Biobase::assayDataElementNames(object)) {
    msg <- c(msg, "The assayData slot must contain a matrix named 'abund'")
  }
  
  if(inherits(object@strategies, "list")) {
    validBcs <- unlist(sapply(names(object), function(bcs) {
      inherits(getStrat(object, bcs), "BiclusterStrategy")
    }))
    if (!all(validBcs)) {
      msg <- c(msg, "All strategies must be BiclusterStrategy objects.")
    } else {
      sapply(names(object), function(bcs) { # Check validity of all strategies
        res <- validObject(getStrat(object, bcs), test = TRUE)
        if(inherits(res, "character")) { msg <<- c(msg, res) }
      })
    }
  } else {
    msg <- c(msg, "The strategies slot must be a 'list' object")
  }
  if(!inherits(object@distance, "dist")) {
    msg <- c(msg, paste("The distance slot must be a 'dist' object as returned",
                        "by dist()"))
  }
  if(is.null(msg)) TRUE else msg
}
setValidity("BiclusterExperiment", validBiclusterExperiment)

#### addStrat ####
#' Add a BiclusterStrategy to a BiclusterExperiment
#' 
#' Due to requirements of various biclustering methods, this function may warn
#' that the chosen biclustering method was overridden.
#' 
#' Also, if any matrix elements of abund(BiclusterExperiment) are missing, then
#' the BiclusterExperiment returned is not guaranteed to have the same number of
#' samples and features.
#' 
#' @export
setGeneric("addStrat", function(bce, k, 
                                method,
                                maxNa = 1) {
  standardGeneric("addStrat")
})

setMethod("addStrat", c(bce = "BiclusterExperiment", k = "numeric", 
                        method = "character"), 
          function(bce, k, method = c("als-nmf", "svd-pca", "snmf",
                                      "nipals-pca", "plaid"), maxNa) {
            # Validate parameters
            if(!is.wholenumber(k)) {
              stop("Arg \"k\" must be a whole number.")
            }
            
            # k must be whole number, smaller than both dimensions of m
            m <- t(as.matrix(bce))
            tryCatch(k <- validateKM(k, m),
                     error = function(e) {
                       warning(paste("Initializing k to the size of the smaller matrix",
                                     "dimension."))
                       k <<- min(nrow(m), ncol(m))
                     }
            )
            
            if(!(maxNa <= 1 && maxNa >= 0)) {
              stop("Arg \"maxNa\" must be in the range of 0 to 1.")
            }
            
            bce <- clean(bce, maxNa)
            silent <- FALSE
            method.orig <- method
            if(length(method) > 1) {
              method <- "als-nmf"
              silent <- TRUE # User does not need any warnings regarding algorithm choice
              # (suppressing warnings is possible only because the BiclusterStrategy
              # constructor does not return any user-informative warnings.
            }
            method <- match.arg(method)
            
            tryCatch({
              mat <- t(as.matrix(bce))
              if(silent) {
                bcs <- suppressWarnings(BiclusterStrategy(m = mat, 
                                                          k = k, method = method))
              } else {
                bcs <- BiclusterStrategy(m = mat, k = k, method = method)
              }
            }, error = function(e) {
              maxNa <- maxNa - (maxNa / 2)
              message(paste("Cleaning with maxNAs at", maxNa))
              
              # Call recursively until success.
              addStrat(bce, k, method.orig, maxNa)
            })
            
            name <- name(bcs)
            bce@strategies[[name]] <- bcs
            if(validObject(bce)) {
              message(paste("Added BiclusterStrategy named", name))
              bce
            }
          })

# FIXME adapt from ExpressionSet method exprs so environment-style assayData can be accessed
#' @export
setMethod("as.matrix", "BiclusterExperiment", function(x) {
  Biobase::assayDataElement(x, "abund")
})

setMethod("clean", c(object = "BiclusterExperiment"), function(object, maxNa) {
  results <- clean(t(as.matrix(object)), maxNa, TRUE)
  # [[2]] contains a vector of indexes of the remaining columns
  # [[1]] contains the cleaned matrix itself
  bce <- object[results[[2]][[2]], results[[2]][[1]]]
  distance(bce) <- dist(results[[1]], method = "euclidean")
  strategies(bce) <- list()
  
  if(validObject(bce, test = FALSE)) bce
})

setMethod("distance", c(bce = "BiclusterExperiment"), function(bce) {
  as.matrix(bce@distance)
}
)
setGeneric("distance<-", function(object, value) standardGeneric("distance<-"))
setReplaceMethod("distance", signature(object = "BiclusterExperiment",
                                       value = "dist"),
                 function(object, value) { object@distance <- value 
                 return(object)}
)

setMethod("getStrat", c(bce = "BiclusterExperiment"), function(bce, id) {
  # if (inherits(id, "character")) { 
  bce@strategies[[id]] 
  # }
  # else {bce@strategies[idL]}
})

#' Names of BiclusterStrategies in this BiclusterExperiment
#' @export
setMethod("names", "BiclusterExperiment", function(x) names(x@strategies))

#' @export
setGeneric("strategies", signature = "bce", function(bce) {
  standardGeneric("strategies")
})
setMethod("strategies", c(bce = "BiclusterExperiment"), function(bce) {
  bce@strategies
})

setGeneric("strategies<-", function(object, value) standardGeneric("strategies<-"))
setReplaceMethod("strategies", signature(object = "BiclusterExperiment",
                                         value = "list"),
                 function(object, value) { 
                   object@strategies <- value 
                   if(validObject(object, test = FALSE)) {
                     object }
                 }
)

#### Abundance heatmap #########################################################
#' Abundance heatmap
#'
#' Plot a heatmap from values in a BiclusterExperiment object.
#'
#' This method implements /code{pheatmap::pheatmap} for BiclusterExperiment objects.
#' 
#' Strategy and cluster must be given if ordering = "cluster" or annotateCluster = TRUE.
#' 
#' @export
setMethod("plot", c(x = "BiclusterExperiment"), 
          function(x, logBase = 0, phenoLabels = c(), biclustLabels = c(), 
                   ordering = c("input", "distance", "cluster"), strategy = "", 
                   rowNames = FALSE, colNames = FALSE) {
            data <- t(as.matrix(x))
            if(logBase != 0) {
              signs <- sign(data)
              data <- log(abs(data), logBase) * signs
              data[is.nan(data)] <- 0
            }
            
            # pheatmap throws cryptic error without rownames
            if(is.null(row.names(data))) row.names(data) <- seq_len(nrow(data))
            
            # Validate requested annotations and parse into a dataframe
            annots <- createAnnots(x, rownames(data), strategy, phenoLabels, biclustLabels)
            
            ordering <- match.arg(ordering)
            
            # FIXME prevent shuffling of annotation colors
            if (ordering == "input") {
              cluster_rows <- FALSE
              cdist <- NULL
              # ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, cluster_cols = FALSE,
              # show_rownames = rowNames, show_colnames = colNames, annotation_row = annots, legend = legend)
            } else if (ordering == "distance") {
              cluster_rows <- TRUE
              cdist <- x@distance
              # ph <- pheatmap::pheatmap(data, clustering_distance_rows = x@distance, cluster_cols = FALSE, show_rownames = rowNames, show_colnames = colNames, annotation_row = annots, legend = legend)
            } else if(ordering == "cluster") {
              cluster_rows <- TRUE
              cdist <- dist(pred(getStrat(x, strategy)), method = "euclidean")
              # ph <- pheatmap::pheatmap(data, cluster_rows = TRUE, 
              # clustering_distance_rows = clusterDist,
              # cluster_cols = FALSE, 
              # show_rownames = rowNames, 
              # show_colnames = colNames, 
              # annotation_row = annots, legend = legend)
            }
            
            ph <- pheatmap::pheatmap(data, cluster_rows = cluster_rows, 
                                     clustering_distance_rows = cdist,
                                     cluster_cols = FALSE,
                                     show_rownames = rowNames, show_colnames = colNames, annotation_row = annots)
          }
)

#### Distance heatmap #########################################################
#' @export
setGeneric("plotDist", signature = "x", function(x, ...) {
  standardGeneric("plotDist")
})

#' Distance heatmap
#'
#' Plot a heatmap from distance values in the BiclusterExperiment object
#'
#' @rdname plotDist
#' @aliases plotDist
setMethod("plotDist", signature(x = "BiclusterExperiment"), 
          function(x, distType = c("euclidean", "pearson"), phenoLabels = c(), biclustLabels = c(), 
                   ordering = c("input", "distance", "cluster"), strategy = "", 
                   rowColNames = FALSE) {
            if(distType == "euclidean") { data <- distance(x) }
            else {data <- 1 - cor(t(as.matrix(x)), method = "pearson")}
            # pheatmap throws cryptic error without rownames
            row.names(data) <- seq_len(nrow(data))
            # Validate requested annotations and parse into a dataframe
            annots <- createAnnots(x, rownames(data), strategy, phenoLabels, biclustLabels)
            
            ordering <- match.arg(ordering)
            
            if (ordering == "input") {
              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, cluster_cols = FALSE,
                                       show_rownames = rowColNames, show_colnames = rowColNames, annotation_row = annots)
            } else if (ordering == "distance") {
              ph <- pheatmap::pheatmap(data, clustering_distance_rows = x@distance, clustering_distance_cols = x@distance, show_rownames = rowColNames, show_colnames = rowColNames, annotation_row = annots)
            } else if(ordering == "cluster") {
              clusterDist <- dist(pred(getStrat(x, strategy)), method = "euclidean")
              ph <- pheatmap::pheatmap(data, clustering_distance_rows = clusterDist,
                                       clustering_distance_cols = clusterDist, 
                                       show_rownames = rowColNames, 
                                       show_colnames = rowColNames, 
                                       annotation_row = annots)
            }
          })

#### Cluster stability plot ####################################################
#' @export
setGeneric("plotClustStab", signature = "obj", function(obj, ...) {
  standardGeneric("plotClustStab")
})

#' Cluster stability plot
#'
#' Plots the stability of clusters in the given BiclusterExperiment.
#'
#' @rdname plotClustStab
#' @aliases plotClustStab
setMethod(
  "plotClustStab", signature(obj = "BiclusterExperiment"),
  function(obj, clusters) {
    if (length(obj@strategies) == 0) {
      warning(paste0("Please calculate clusters first!"))
      return(obj)
    }
    
    # calculate consensus too?
    
    # calculate stability of the clusters check if there are more than 1 k value in ks range
    stability <- NULL
    # stability <- calculate_stability(metadata(object)$sc3$consensus, k)
    stability <- abs(runif(n = clusters, min = 0, max = 1))
    
    d <- data.frame(Cluster = factor(1:length(stability)), Stability = stability)
    ggplot2::ggplot(d, ggplot2::aes_string(x = "Cluster", y = "Stability")) +
      ggplot2::geom_bar(stat = "identity") + ggplot2::ylim(0, 1) +
      ggplot2::labs(x = "Cluster", y = "Stability Index") + ggplot2::theme_bw()
  }
)

#### Factor matrix heatmap ###################################################
#' @export
setGeneric("heatmapFactor", signature = "obj", function(obj, ...) {
  standardGeneric("heatmapFactor")
})

#' Factor matrix heatmap
#' 
#' Display scores for all clusters in one heatmap
setMethod("heatmapFactor", c(obj = "BiclusterExperiment"), 
          function(obj, strategy = "", type = c("score", "loading"), phenoLabels = c(), biclustLabels = c(), 
                   ordering = c("input", "distance", "cluster"), 
                   colNames = FALSE) {
            type <- match.arg(type)
            if(type == "score") {
              data <- t(score(getStrat(obj, strategy)))
              annots <- createAnnots(obj, colnames(data), strategy, phenoLabels, biclustLabels)
            } else {
              data <- loading(getStrat(obj, strategy))
              annots <- NA
            }
            
            if (is.null(rownames(data))) {
              rownames(data) <- sapply(seq_len(nrow(data)), function(x) {
                paste0("Bicluster.", x)
              }
              )
            }
            if(is.null(colnames(data))) {
              colnames(data) <- seq_len(ncol(data))
            }
            # Validate requested annotations and parse into a dataframe
            
            ordering <- match.arg(ordering)
            silent = TRUE
            if (ordering == "input") {
              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, 
                                       cluster_cols = FALSE, 
                                       show_colnames = colNames, 
                                       annotation_col = annots, silent = silent)
            } else if (ordering == "distance") {
              distance <- dist(as.matrix(obj), method = "euclidean")
              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE,
                                       clustering_distance_cols = distance,
                                       show_colnames = colNames,
                                       annotation_col = annots, silent = silent)
            } else {
              # FIXME this is not valid for loadings
              clusterDist <- dist(pred(getStrat(obj, strategy)), method = "euclidean")
              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, 
                                       clustering_distance_cols = clusterDist,
                                       show_colnames = colNames, 
                                       annotation_col = annots, silent = silent)
            }
            ph
          }
)

#### Score plot ################################################################
#' @param thresholds A single numeric, vector, or matrix of thresholds. See 
#'   documentation
#' @export
setGeneric("plotSamples", signature = "obj", function(obj, thresholds = NULL, strategy, bicluster, ordering) {
  standardGeneric("plotSamples")
})


#' Score plot
#'
#' Plot cluster membership scores
#'
#' @rdname plotSamples
#' @aliases plotSamples
setMethod(
  "plotSamples", signature(obj = "BiclusterExperiment"),
  function(obj, thresholds, strategy = "", bicluster = "Bicluster.1", 
           ordering = c("input", "distance", "cluster")) {
    bcs <- getStrat(obj, strategy)
    bicluster <- validateBiclustNames(biclustNames = bicluster, bcs = bcs)
    if(length(bicluster) != 1) {
      stop("Argument \"bicluster\" was not valid.")
    }
    # User is allowed to provide custom threshold when this function is called outside GUI
    if (is.null(thresholds)) {
      thresholds <- bcs@scoreThresh[bicluster, ]
      names(thresholds) <- colnames(bcs@scoreThresh)
      # FIXME Try changing the BCE/BCS constructors so this check isn't necessary
    } else if (inherits(thresholds, "numeric")) {
      names(thresholds) <- as.character(thresholds)
    } else {
      stop("Argument \"thresholds\" must be numeric. Use this argument only if
           you want to replce mfBiclust thresholds with your own.")
    }
    
    # Define colors
    cols <- RColorBrewer::brewer.pal(
      if (length(thresholds) < 3) {3}
      else {length(thresholds)}, 
      "Dark2"
    )
    
    data <- score(bcs)[, bicluster]
    ordering <- match.arg(ordering)
    # Plot
    if(ordering == "input") {
      plot(1:nrow(t(as.matrix(obj))), data,
           xlab = "Sample", ylab = "Score", main = bicluster)
    } else if (ordering == "distance") {
      ord <- hclust(obj@distance)$order
      plot(1:nrow(t(as.matrix(obj))), data[ord],
           xlab = "Sample", ylab = "Score", xaxt = "n", main= bicluster)
      axis(1, at = 1:nrow(t(as.matrix(obj))), labels = as.character(ord))
    } else if (ordering == "cluster") {
      ord <- order(data, decreasing = TRUE)
      plot(1:nrow(t(as.matrix(obj))), data[ord],
           xlab = "Sample", ylab = "Score", xaxt = "n", main= bicluster)
      axis(1, at = 1:nrow(t(as.matrix(obj))), labels = as.character(ord))
    }
    
    mapply(function(y, color) {abline(h = y, col = color, lwd = 2)},
           y = thresholds, color = cols[seq_along(thresholds)]
    )
    
    legend("topright", legend = capitalize(names(thresholds)), col = cols[seq_along(names(thresholds))], lty = 1, 
           lwd = 2, cex = 0.8
    )
  }
)

#### Biomarker plot ############################################################
#' Biomarker plot
#'
#' Plot features with some measure of relevance on the y-axis
#'
#' @rdname plotMarkers
#' @aliases plotMarkers
#' @export
setGeneric("plotMarkers", signature = "obj", function(obj, thresholds = NULL, ...) {
  standardGeneric("plotMarkers")
})

#' cluster must be numeric
#' @export
setMethod("plotMarkers", signature(obj = "BiclusterExperiment"),
          function(obj, thresholds, strategy = "", bicluster = "Bicluster.1", 
                   ordering = c("input", "cluster")) {
            
            bcs <- getStrat(obj, strategy)
            bicluster <- validateBiclustNames(biclustNames = bicluster, bcs = bcs)
            if(length(bicluster) != 1) {
              stop("Argument \"bicluster\" was not valid.")
            }
            # User is allowed to provide custom threshold when this function is called outside GUI
            if (is.null(thresholds)) {
              thresholds <- bcs@loadingThresh[bicluster, ]
              names(thresholds) <- colnames(bcs@loadingThresh)
            } else if (inherits(thresholds, "numeric")) {
              names(thresholds) <- as.character(thresholds)
            } else {
              stop("Argument \"thresholds\" must be numeric. Use this argument only if
                   you want to replce mfBiclust thresholds with your own.")
            }
            
            # Define colors
            cols <- RColorBrewer::brewer.pal(
              if (length(thresholds) < 3) {3}
              else {length(thresholds)}, 
              "Dark2"
            )
            
            data <- t(loading(bcs))[, bicluster]
            ordering <- match.arg(ordering)
            
            # FIXME ADD ACCURACY, RECOVERY, etc. on
            # PLOT. ALLOW TO HIGHLIGHT KNOWN FEATURES.
            # Plot
            if(ordering == "input") {
              plot(1:ncol(t(as.matrix(obj))), data,
                   xlab = "Feature", ylab = "Loading", main = bicluster)
            } else if (ordering == "distance") {
              ord <- hclust(dist(t(as.matrix(obj)), method = "euclidean"))$order
              plot(1:ncol(t(as.matrix(obj))), data[ord],
                   xlab = "Feature", ylab = "Loading", xaxt = "n", main= bicluster)
              axis(1, at = 1:ncol(t(as.matrix(obj))), labels = as.character(ord))
            } else if (ordering == "cluster") {
              ord <- order(data, decreasing = TRUE)
              plot(1:ncol(t(as.matrix(obj))), data[ord],
                   xlab = "Feature", ylab = "Loading", xaxt = "n", main= bicluster)
              axis(1, at = 1:ncol(as.matrix(obj)), labels = as.character(ord))
            }
            
            mapply(function(y, color) {abline(h = y, col = color, lwd = 2)},
                   y = thresholds, color = cols[seq_along(thresholds)]
            )
            
            legend("topright", legend = capitalize(names(thresholds)), col = cols[seq_along(names(thresholds))], lty = 1, 
                   lwd = 2, cex = 0.8
            )
          }
)

#### PCA Plot ##################################################################
#' Biomarker plot
#'
#' Plot features with some measure of relevance on the y-axis
#'
#' @rdname plotMarkers
#' @aliases plotMarkers
#' @export
setGeneric("pca", signature = "bce", function(bce) {
  standardGeneric("pca")
})
setMethod("pca", signature(bce = "BiclusterExperiment"), function(bce) {
  pcaStrats <- unlist(lapply(strategies(bce), function(strat) {
    (method(strat) == "nipals-pca" || method(strat) == "svd-pca")
  }))
  if(any(pcaStrats)) {
    # PCA has already been performed
    strat <- getStrat(bce, min(which(pcaStrats))) # A PCA-encapsulating strategy
    var <- nrow(as.matrix(bce))
    xs <- strat@factors@fit@W[, 1]
    var1 <- round(sd(xs) ^ 2 / var, 1)
    ys <- strat@factors@fit@W[, 2]
    var2 <- round(sd(ys) ^ 2 / var, 1)
  } else {
    m <- as.matrix(bce)
    prcmp <- prcomp(m)
    xs <- prcmp$x[, 1]
    var1 <- round(prcmp$sdev[1]^2 / sum(prcmp$sdev^2) * 100, 1)
    ys <- prcmp$x[, 2]
    var2 <- round(prcmp$sdev[2]^2 / sum(prcmp$sdev^2) * 100, 1)
  }
  plot(xs, ys, pch = 16,
       xlab = paste0("PC1: ", 
                     var1, 
                     "% variance"),
       ylab = paste0("PC2: ", 
                     var2, 
                     "% variance"))
})

