#' @include BiclusterExperiment.R
#' @include BiclusterStrategy.R
#' @include helperFunctions.R

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
setGeneric("plotDist", signature = "bce", function(bce, type, distType, ordering,
                                                ...) {
  standardGeneric("plotDist")
})

#' Distance heatmap
#'
#' Plot a heatmap from distance values in the BiclusterExperiment object
#'
#' @rdname plotDist
#' @aliases plotDist
setMethod(
  "plotDist", signature(bce = "BiclusterExperiment"), 
  function(bce, type = c("samples", "features"),
           distType = c("euclidean", "pearson"), 
           ordering = c("input", "distance", "cluster"), phenoLabels = c(),
           biclustLabels = c(), strategy = "", rowColNames = FALSE) {
    type <- match.arg(type)
    ordering <- match.arg(ordering)
    
    m <- as.matrix(bce) # get data
    m <- if((distType == "euclidean") == (distType == "samples")) {
      m 
    } else {
      t(m) # transpose if necessary
    }
    if(distType == "euclidean") {
      dObj <- dist(m)
      data <- as.matrix(dObj)
    } else if(distType == "pearson") {
      data <- 1 - cor(m)
      # Sometimes this results in redundant t() calls
      if(ordering == "distance") dObj <- dist(t(m))
    }
    
    # pheatmap throws cryptic error without rownames
    row.names(data) <- seq_len(nrow(data))
    # Validate requested annotations and parse into a dataframe
    annots <- createAnnots(bce, rownames(data), strategy, phenoLabels, biclustLabels)
    
    if (ordering == "input") {
      ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, cluster_cols = FALSE,
                               show_rownames = rowColNames, show_colnames = rowColNames, annotation_row = annots)
    } else if (ordering == "distance") {
      clust <- hclust(dObj)
      ph <- pheatmap::pheatmap(data, cluster_rows = clust,
                               cluster_cols = clust,
                               show_rownames = rowColNames,
                               show_colnames = rowColNames,
                               annotation_row = annots)
    } else if(ordering == "cluster") {
      clusFunction <- if(type == "samples") {
        clusteredSamples 
      } else {
        clusteredFeatures
      }
      clus <- hclus(dist(clusFunction(getStrat(bce, strategy)),
                          method = "euclidean"))
      ph <- pheatmap::pheatmap(data, cluster_rows = clus,
                               cluster_cols = clus, 
                               show_rownames = rowColNames, 
                               show_colnames = rowColNames, 
                               annotation_row = annots)
    }
  })
#' 
#' #### Cluster stability plot ####################################################
#' #' @export
#' setGeneric("plotClustStab", signature = "obj", function(obj, ...) {
#'   standardGeneric("plotClustStab")
#' })
#' 
#' #' Cluster stability plot
#' #'
#' #' Plots the stability of clusters in the given BiclusterExperiment.
#' #'
#' #' @rdname plotClustStab
#' #' @aliases plotClustStab
#' setMethod(
#'   "plotClustStab", signature(obj = "BiclusterExperiment"),
#'   function(obj, clusters) {
#'     if (length(obj@strategies) == 0) {
#'       warning(paste0("Please calculate clusters first!"))
#'       return(obj)
#'     }
#'     
#'     # calculate consensus too?
#'     
#'     # calculate stability of the clusters check if there are more than 1 k value in ks range
#'     stability <- NULL
#'     # stability <- calculate_stability(metadata(object)$sc3$consensus, k)
#'     stability <- abs(runif(n = clusters, min = 0, max = 1))
#'     
#'     d <- data.frame(Cluster = factor(1:length(stability)), Stability = stability)
#'     ggplot2::ggplot(d, ggplot2::aes_string(x = "Cluster", y = "Stability")) +
#'       ggplot2::geom_bar(stat = "identity") + ggplot2::ylim(0, 1) +
#'       ggplot2::labs(x = "Cluster", y = "Stability Index") + ggplot2::theme_bw()
#'   }
#' )

#### Factor matrix heatmap ###################################################
#' @export
setGeneric("factorHeatmap", signature = c("bce", "bcs"), function(bce, bcs,
                                                                  type, ...) {
  standardGeneric("factorHeatmap")
})
setMethod("factorHeatmap", c(bce = "BiclusterExperiment", bcs = "character"),
          function(bce, bcs, type, ...) {
            factorHeatmap(bce, getStrat(bce, bcs), type, ...)
          })
setMethod("factorHeatmap", c(bce = "BiclusterExperiment", bcs = "numeric"),
          function(bce, bcs, type, ...) {
            factorHeatmap(bce, getStrat(bce, bcs), type, ...)
          })
#' Factor matrix heatmap
#' 
#' Display scores for all clusters in one heatmap
setMethod(
  "factorHeatmap", c(bce = "BiclusterExperiment", bcs = "BiclusterStrategy"),
  function(
    bce, bcs, type = c("score", "loading"), phenoLabels = c(),
    biclustLabels = c(), ordering = c("input", "distance", "cluster"),
    colNames = FALSE) {
    type <- match.arg(type)
    
    if(type == "score") {
      data <- t(score(bcs))
      # Validate requested annotations and parse into a dataframe
      annots <- createAnnots(bce, colnames(data), bcs, phenoLabels, biclustLabels)
    } else {
      data <- loading(bcs)
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
    ordering <- match.arg(ordering)
    silent = TRUE
    if (ordering == "input") {
      ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, 
                               cluster_cols = FALSE, 
                               show_colnames = colNames, 
                               annotation_col = annots, silent = silent)
    } else if (ordering == "distance") {
      distance <- dist(as.matrix(bce), method = "euclidean")
      ph <- pheatmap::pheatmap(data, cluster_rows = FALSE,
                               clustering_distance_cols = distance,
                               show_colnames = colNames,
                               annotation_col = annots, silent = silent)
    } else { # cluster reordering
      clusFunction <- if(type == "score") {
        clusteredSamples
      } else {
        clusteredFeatures
      }
      clusterDist <- dist(clusFunction(bcs), method = "euclidean")
      ph <- pheatmap::pheatmap(data, cluster_rows = FALSE, 
                               clustering_distance_cols = clusterDist,
                               show_colnames = colNames, 
                               annotation_col = annots, silent = silent)
    }
    ph
  }
)

#### Threshold plot ################################################################
#' @param thresholds A single numeric, vector, or matrix of thresholds. See 
#'   documentation
#' @export
setGeneric("plotThreshold", signature = c("bce", "bcs"),
           function(bce, bcs, type, bicluster, ...) {
             standardGeneric("plotThreshold")
           })
setMethod("plotThreshold",
          signature(bce = "BiclusterExperiment", bcs = "character"),
          function(bce, bcs, type, bicluster, ...) {
            plotThreshold(bce = bce, bcs = getStrat(bce, bcs), type, bicluster, ...)
          }
)
setMethod("plotThreshold",
          signature(bce = "BiclusterExperiment", bcs = "numeric"),
          function(bce, bcs, type, bicluster, ...) {
            plotThreshold(bce, getStrat(bce, bcs), type, bicluster, ...)
          }
)
#' Threshold plot
#'
#' Plot cluster membership values and thresholding
#'
#' @rdname plotThreshold
setMethod(
  "plotThreshold", signature(bce = "BiclusterExperiment", bcs = "BiclusterStrategy"),
  function(bce, bcs, type = c("score", "loading"), bicluster = "Bicluster.1", 
           thresholds = NULL, ordering = c("input", "distance", "cluster"),
           xlabs = FALSE) {
    bicluster <- validateBiclustNames(biclustNames = bicluster, bcs = bcs)

    if(length(bicluster) != 1) {
      stop("Argument \"bicluster\" was not valid.")
    }
    # User is allowed to provide custom threshold when this function is called outside GUI
    if (is.null(thresholds)) {
      thresholds <- get(paste0(type, "Thresh"))(bcs)[bicluster]
      names(thresholds) <- threshAlgo(bcs)
    } else if (inherits(thresholds, "numeric")) {
      names(thresholds) <- as.character(thresholds)
    } else {
      stop(paste("Argument \"thresholds\" must be numeric or NULL Use this",
                 "argument only if you want to replce mfBiclust thresholds",
                 "with your own."))
    }
    
    # Define colors
    cols <- RColorBrewer::brewer.pal(
      if (length(thresholds) < 3) {3}
      else {length(thresholds)}, 
      "Dark2"
    )
    data <- if(type == "score") { get(type)(bcs) } else {
      t(get(type)(bcs))
    }
    data <- data[, bicluster, drop = TRUE]
    
    ordering <- match.arg(ordering)
    # Plot
    xs <- seq_along(data)
    xlab <- if(type == "score") "Samples" else "Features"
    ylab = capitalize(type)
    
    if(ordering == "input") {
      ord <- xs
    } else if (ordering == "distance") {
      m <- if(type == "score") t(as.matrix(bce)) else as.matrix(bce)
      ord <- hclust(dist(m))$order
    } else if (ordering == "cluster") {
      clusFunction <- if(type == "score") {
        function(bcs) clusteredSamples(bcs)
      } else {
        function(bcs) t(clusteredFeatures(bcs))
      }
      ord <- hclust(dist(clusFunction(bcs)))$order
    }
    plot(xs, data[ord], xlab = xlab, ylab = ylab, xaxt = "n", main = bicluster)
    axis(1, at = xs, labels = if(xlabs) {
      names(data[ord])
      } else {
        rep("", length(xs))
        }
      )
    
    mapply(function(y, color) {abline(h = y, col = color, lwd = 2)},
           y = thresholds, color = cols[seq_along(thresholds)]
    )
    
    legend("topright", legend = capitalize(names(thresholds)), col = cols[seq_along(names(thresholds))], lty = 1, 
           lwd = 2, cex = 0.8
    )
    }
)

#### PCA Plot ##################################################################
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
