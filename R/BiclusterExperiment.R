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

#### CONSTRUCTORS #############################################################
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
#' @return an instance of BiclusterExperiment-class
#'   containing the following slots, accessed with @@: 
#'   data: Object of class \code{\link{matrix}}. The original data. 
#'   annot: Object of class \code{\link{data.frame}}. Annotations provided by the user
#'   strategies: A \code{\link{list}} of \code{BiclusterStrategy} objects
#'   distance: Object of class \code{\link{matrix}}. Pairwise distances
#'     between the rows of \code{data}
setGeneric("BiclusterExperiment", function(bcs, ...) {
  standardGeneric("BiclusterExperiment")
})

#' Careful, for m, rows are samples and columns are features
#' eSet objects store assayData transposed: rows are features and columns are samples.
#' For this reason I wrote a getter that returns a matrix with rows as samples, columns as features.
setMethod("BiclusterExperiment", c(bcs = "BiclusterStrategy"), function(bcs, m = matrix(), annot = data.frame(), bcv = FALSE) {
  if (bcv == TRUE) {
    warning("Bi-cross-validation is still under development. msNMF
                      cannot predict the optimal clustering strategy.")
  }
  
  d <- as.matrix(dist(m, method = "euclidean"))
  
  if(is.null(row.names(annot)) && is.null(row.names(m))) {
    # If no sample names provided, just call them Sample.1, Sample.2, Sample.3
    row.names(m) <- unlist(sapply(seq_along(nrow(m)), function(s) {
      paste0("Sample.", s)
    }
    ))
    row.names(annot) <- row.names(m)
  } else if(is.null(row.names(annot))) {
    # If sample names provided in one argument but not the other, just assume
    row.names(annot) <- row.names(m)
  } else if(is.null(row.names(m))) {
    row.names(m) <- row.names(annot)
  }
  
  ad <- Biobase::assayDataNew(storage.mode = "list")
  ad[[1]] <- t(m)
  
  bcs <- list(bcs)
  names(bcs) <- name(bcs[[1]])
  
  new("BiclusterExperiment", assayData = t(m), phenoData = AnnotatedDataFrame(data = annot), strategies = bcs, distance = d)
})

setMethod("BiclusterExperiment", c(bcs = "list"), function(bcs, m = matrix(), annot = data.frame(), bcv = FALSE) {
  if (!all(unlist(lapply(bcs, function(obj) inherits(obj, "BiclusterStrategy"))))) {
    stop("Argument \"bcs\" must contain only BiclusterStrategy objects.")
  }
  if (length(bcs) == 0) {
    warning("Since argument \"bcs\" is an empty list, an empty
                    BiclusterExperiment object will be instantiated.")
  }
  if (bcv == TRUE) {
    warning("Bi-cross-validation is still under development. msNMF cannot 
            predict the optimal clustering strategy.")
  }
  
  d <- dist(m, method = "euclidean")
  
  if(is.null(row.names(annot)) && is.null(row.names(m))) {
    # If no sample names provided, just call them Sample.1, Sample.2, Sample.3
    row.names(m) <- unlist(sapply(seq_along(nrow(m)), function(s) {
      paste0("Sample.", s)
    }
    ))
    row.names(annot) <- row.names(m)
  } else if(is.null(row.names(annot))) {
    # If sample names provided in one argument but not the other, just assume
    row.names(annot) <- row.names(m)
  } else if(is.null(row.names(m))) {
    row.names(m) <- row.names(annot)
  }
  
  ad <- Biobase::assayDataNew(storage.mode = "list")
  ad[[1]] <- t(m)
  
  names(bcs) <- lapply(bcs, function(bcs) {name(bcs)})
  
  new("BiclusterExperiment", assayData = ad, phenoData = AnnotatedDataFrame(data = annot), strategies = bcs, distance = d)
})

#### METHODS ###################################################################

#' @export
setMethod("as.matrix", "BiclusterExperiment", function(x) t(assayData(x)[[1]]))

#' Names of BiclusterStrategies in this BiclusterExperiment
#' @export
setMethod("names", "BiclusterExperiment", function(x) names(x@strategies))

setMethod("getStrat", c(bce = "BiclusterExperiment"), function(bce, stratName) {
  bce@strategies[[stratName]]
}
)

setMethod("distMat", c(bce = "BiclusterExperiment"), function(bce) {
  as.matrix(bce@distance)
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
            
            data <- as.matrix(x)
            if(logBase != 0) {
              signs <- sign(data)
              data <- log(abs(data), logBase) * signs
              data[is.nan(data)] <- 0
            }
            # pheatmap throws cryptic error without rownames
            row.names(data) <- seq_len(nrow(data))
            # Validate requested annotations and parse into a dataframe
            annots <- createAnnots(x, rownames(data), strategy, phenoLabels, biclustLabels)
            
            ordering <- match.arg(ordering)
            
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

# adapted from sc3
# make_ann_df_for_heatmaps <- function(annotations, showAnnot) {
#   if (any(!showAnnot %in% colnames(annotations))) {
#     # feedback; irrelevant for GUI
#     warning(paste0(
#       "Provided annotation tracks '",
#       paste(showAnnot[!showAnnot %in% colnames(annotations)], 
#             collapse = "', '"),
#       "' do not exist in the dataframe. They will not be displayed."
#     ))
#     showAnnot <- showAnnot[showAnnot %in% colnames(annotations)]
#   }
#   ann <- annotations[showAnnot]
#   
#   # remove columns with 1 value only
#   if (length(showAnnot) > 1) {
#     if (ncol(ann) == 0) {
#       ann <- NULL
#     } else {
#       ann <- as.data.frame(lapply(ann, function(x) {
#         if (nlevels(as.factor(x)) > 9) {
#           x
#         } else {
#           as.factor(x)
#         }
#       }))
#     }
#   } else {
#     ann <- as.data.frame(ann)
#     colnames(ann) <- showAnnot
#     # ann <- as.data.frame(lapply(ann, function(x) {
#     #   if (nlevels(as.factor(x)) > 9) {
#     #     return(x)
#     #   } else {
#     #     return(as.factor(x))
#     #   }
#     # }))
#   }
#   return(ann)
# }

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
            if(distType == "euclidean") { data <- distMat(x) }
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
            if (ordering == "input") {
              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, 
                                       cluster_cols = FALSE, 
                                       show_colnames = colNames, 
                                       annotation_col = annots)
            } else if (ordering == "distance") {
              distance <- dist(as.matrix(obj), method = "euclidean")
              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE,
                                       clustering_distance_cols = distance,
                                       show_colnames = colNames,
                                       annotation_col = annots)
            } else {
              # FIXME this is not valid for loadings
              clusterDist <- dist(pred(getStrat(obj, strategy)), method = "euclidean")
              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, 
                                       clustering_distance_cols = clusterDist,
                                       show_colnames = colNames, 
                                       annotation_col = annots)
            }
}
)

#### Score plot ################################################################
#' @export
setGeneric("plotSamples", signature = "obj", function(obj, thresholds = NULL, ...) {
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
      plot(1:nrow(as.matrix(obj)), data,
           xlab = "Sample", ylab = "Score", main = bicluster)
    } else if (ordering == "distance") {
      ord <- hclust(obj@distance)$order
      plot(1:nrow(as.matrix(obj)), data[ord],
           xlab = "Sample", ylab = "Score", xaxt = "n", main= bicluster)
      axis(1, at = 1:nrow(as.matrix(obj)), labels = as.character(ord))
    } else if (ordering == "cluster") {
      ord <- order(data, decreasing = TRUE)
      plot(1:nrow(as.matrix(obj)), data[ord],
           xlab = "Sample", ylab = "Score", xaxt = "n", main= bicluster)
      axis(1, at = 1:nrow(as.matrix(obj)), labels = as.character(ord))
    }

    mapply(function(y, color) {abline(h = y, col = color, lwd = 2)},
           y = thresholds, color = cols[seq_along(thresholds)]
    )
    
    legend("topright", legend = names(thresholds), col = cols[seq_along(names(thresholds))], lty = 1, 
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
          
            # Plot
            if(ordering == "input") {
              plot(1:ncol(as.matrix(obj)), data,
                   xlab = "Feature", ylab = "Loading", main = bicluster)
            } else if (ordering == "distance") {
              ord <- hclust(dist(t(as.matrix(obj)), method = "euclidean"))$order
              plot(1:ncol(as.matrix(obj)), data[ord],
                   xlab = "Feature", ylab = "Loading", xaxt = "n", main= bicluster)
              axis(1, at = 1:ncol(as.matrix(obj)), labels = as.character(ord))
            } else if (ordering == "cluster") {
              ord <- order(data, decreasing = TRUE)
              plot(1:ncol(as.matrix(obj)), data[ord],
                   xlab = "Feature", ylab = "Loading", xaxt = "n", main= bicluster)
              axis(1, at = 1:ncol(as.matrix(obj)), labels = as.character(ord))
            }
            
            mapply(function(y, color) {abline(h = y, col = color, lwd = 2)},
                   y = thresholds, color = cols[seq_along(thresholds)]
            )
            
            legend("topright", legend = names(thresholds), col = cols[seq_along(names(thresholds))], lty = 1, 
                   lwd = 2, cex = 0.8
            )
            }
)