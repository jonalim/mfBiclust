#' @include BiclusterExperiment.R
#' @include BiclusterStrategy.R
#' @include helperFunctions.R

#### Abundance heatmap #########################################################
#' @importFrom pheatmap pheatmap
#' @rdname plot
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
              if(is.null(row.names(data))) {
                  row.names(data) <- seq_len(nrow(data))
              }

              # Validate requested annotations and parse into a dataframe
              annots <- createAnnots(x, type = "sample", strategy, phenoLabels,
                                     biclustLabels)

              ordering <- match.arg(ordering)

              # FIXME prevent shuffling of annotation colors?
              if (ordering == "input") {
                  cluster_cols <- FALSE
                  cdist <- NULL
              } else if (ordering == "distance") {
                  cluster_cols <- TRUE
                  cdist <- dist(t(as.matrix(x)))
              } else if(ordering == "cluster") {
                  cluster_cols <- TRUE
                  cdist <- dist(t(clusteredSamples(getStrat(x, strategy))),
                                method = "euclidean")
              }

              ph <- pheatmap::pheatmap(data, cluster_rows = FALSE,
                                       clustering_distance_cols = cdist,
                                       cluster_cols = cluster_cols,
                                       show_rownames = rowNames,
                                       show_colnames = colNames,
                                       annotation_col = annots)
          }
)

#### Distance heatmap #########################################################
#' @return a \code{\link[pheatmap]{pheatmap}-class} object
#' @export
setGeneric("plotDist", signature = "bce", function(
    bce, type, ordering = c("input", "distance", "cluster"),
    distType = c("euclidean", "pearson"), ...) {
    standardGeneric("plotDist")
})

#' Distance heatmap
#'
#' Plot a heatmap of sample distance or feature distance calculated from the raw
#' data in a \code{\link{BiclusterExperiment-class}} object. Default is
#' Euclidean distance between samples.
#'
#' @param bce A \code{BiclusterExperiment} object
#' @param type "sample" or "feature"
#' @param ordering The default \code{ordering = "input"} preserves the order
#'   of samples or features from \code{bce@@assayData}. \code{"distance"}
#'   reorders based on the distance calculated from \code{distType}.
#'   \code{"cluster"} reorders based on bicluster membership.
#' @param distType The formula to use for distance
#' @param ... Optional parameters \code{phenoLabels}, \code{biclustLabels},
#'   \code{strategy}, and \code{rowColNames}, described below:
#' @param phenoLabels an optional character vector of labels to annotate. If
#'   \code{type = "feature"}, \code{phenoLabels} should be column names of
#'   \code{Biobase::phenoData(bce)}
#' @param biclustLabels an optional character vector of labels to annotate.
#'   Should be elements of \code{bcNames(bcs)}. Both \code{phenoLabels} and
#'   \code{biclustLabels} may be specified.
#' @param strategy The name or index of an object in \code{bce@@strategies} from
#'   which to obtain bicluster annotations, or a
#'   \code{\link{BiclusterStrategy-class}} object.
#' @param rowColNames Should feature/sample names be printed
#'
#' @examples
#' bce <- BiclusterExperiment(cancer_benchmark[[1]]$data)
#' plotDist(bce, type = "sample", ordering = "distance", rowColNames = TRUE)
#' plotDist(bce, type = "feature", distType = "pearson", ordering = "distance")
#'
#' @importFrom pheatmap pheatmap
#' @rdname plotDist
#' @aliases plotDist
setMethod(
    "plotDist", signature(bce = "BiclusterExperiment"),
    function(bce, type = c("sample", "feature"),
             ordering = c("input", "distance", "cluster"),
             distType = c("euclidean", "pearson"),
             phenoLabels = c(),
             biclustLabels = c(), strategy = "", rowColNames = FALSE) {
        type <- match.arg(type)
        ordering <- match.arg(ordering)
        distType <- match.arg(distType)

        m <- as.matrix(bce) # get data
        m <- if((distType == "euclidean") == (distType == "sample")) {
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
        annots <- createAnnots(bce, type = "score", strategy, phenoLabels,
                               biclustLabels)

        if (ordering == "input") {
            ph <- pheatmap::pheatmap(data, cluster_rows = FALSE,
                                     cluster_cols = FALSE,
                                     show_rownames = rowColNames,
                                     show_colnames = rowColNames,
                                     annotation_row = annots)
        } else if (ordering == "distance") {
            clust <- hclust(dObj)
            ph <- pheatmap::pheatmap(data, cluster_rows = clust,
                                     cluster_cols = clust,
                                     show_rownames = rowColNames,
                                     show_colnames = rowColNames,
                                     annotation_row = annots)
        } else if(ordering == "cluster") {
            clusFunction <- if(type == "sample") {
                clusteredSamples
            } else {
                clusteredFeatures
            }
            clus <- hclust(dist(clusFunction(getStrat(bce, strategy)),
                                method = "euclidean"))
            ph <- pheatmap::pheatmap(data, cluster_rows = clus,
                                     cluster_cols = clus,
                                     show_rownames = rowColNames,
                                     show_colnames = rowColNames,
                                     annotation_row = annots)
        }
    })

#### Factor matrix heatmap ###################################################
#' Plot a heatmap showing bicluster membership of samples or features
#'
#' Reads data from \code{BiclusterStrategy@@factors} to create a heatmap of
#' bicluster membership across all samples or features.
#'
#' @param bce A BiclusterExperiment object
#' @param bcs The name or index of a BiclusterStrategy contained by \code{bce},
#'   or the BiclusterStrategy object itself
#' @param type either "feature" for feature-bicluster membership or "sample" for
#'   sample-bicluster membership
#' @param ... Optional arguments \code{phenoLabels}, \code{biclustLabels},
#'   \code{ordering}, and \code{colNames}, described below:
#' @param phenoLabels an optional character vector of labels to annotate. If
#'   \code{type = "feature"}, \code{phenoLabels} should be column names of
#'   \code{Biobase::phenoData(bce)}
#' @param biclustLabels an optional character vector of labels to annotate.
#'   Should be elements of \code{bcNames(bcs)}. Both \code{phenoLabels} and
#'   \code{biclustLabels} may be specified.
#' @param ordering The default \code{ordering = "input"} preserves the order
#'   of samples or features from \code{bce@@assayData}. \code{"distance"}
#'   reorders based on Euclidean distance calculated from \code{bce@@assayData}.
#'   \code{"cluster"} reorders based on bicluster membership.
#' @param colNames if \code{TRUE}, labels the samples/features
#'
#' @return a \code{\link[pheatmap]{pheatmap}-class} object
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' factorHeatmap(bce, bcs, type = "sample")
#'
#' @importFrom pheatmap pheatmap
#' @export
setGeneric("factorHeatmap", signature = c("bce", "bcs"),
           function(bce, bcs, type, ordering = "input", ...)
               standardGeneric("factorHeatmap")
)
#' @describeIn factorHeatmap Plots a matrix factor from the
#'   \code{\link{BiclusterStrategy-class}} object named \code{bcs} in
#'   \code{bce@@strategies}.
setMethod("factorHeatmap", c(bce = "BiclusterExperiment", bcs = "character"),
          function(bce, bcs, type, ...) {
              factorHeatmap(bce, getStrat(bce, bcs), type, ordering, ...)
          })
#' @describeIn factorHeatmap Plots a matrix factor from the
#'   \code{\link{BiclusterStrategy-class}} object at the index specified by
#'   \code{bcs}
setMethod("factorHeatmap", c(bce = "BiclusterExperiment", bcs = "numeric"),
          function(bce, bcs, type, ...) {
              factorHeatmap(bce, getStrat(bce, bcs), type, ordering, ...)
          })
#' @describeIn factorHeatmap Plots a matrix factor from \code{bcs}.
setMethod(
    "factorHeatmap", c(bce = "BiclusterExperiment", bcs = "BiclusterStrategy"),
    function(bce, bcs, type = c("feature", "sample"),
             ordering = c("input", "distance", "cluster"), phenoLabels = c(),
             biclustLabels = c(), colNames = FALSE) {

        grDevices::dev.hold()
        on.exit(grDevices::dev.flush())
        plot.new()

        type <- match.arg(type)

        if(type == "feature") {
            data <- t(logicalMatrix2Numeric(featureFactor(bcs)))
            # Validate requested annotations and parse into a dataframe
            annots <- createAnnots(bce, type = type, bcs, phenoLabels, biclustLabels)
        } else {
            data <- logicalMatrix2Numeric(sampleFactor(bcs))
            annots <- createAnnots(bce, type = type, bcs, phenoLabels, biclustLabels)
        }

        if (is.null(rownames(data))) {
            stop("Please ensure all biclusters in the provided BiclusterStrategy are named.")
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
            distance <- dist(t(as.matrix(bce)), method = "euclidean")
            ph <- pheatmap::pheatmap(data, cluster_rows = FALSE,
                                     clustering_distance_cols = distance,
                                     show_colnames = colNames,
                                     annotation_col = annots, silent = silent)
        } else { # cluster reordering
            clusFunction <- if(type == "sample") {
                function(x) t(clusteredSamples(x))
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

#' Quickly preview a matrix as a heatmap
#'
#' Uses R's \code{image()} function to visualize a matrix.
#'
#' @param m A numeric matrix
#'
#' @return Invisibly returns NULL
#'
#' @examples
#' matrixHeatmap(yeast_benchmark[[1]])
#'
#' @export
matrixHeatmap <- function(m) {
    old.par <- par(no.readonly = TRUE)
    par(mar = c(0, 0, 0, 0))
    image(t(apply(m, 2, rev)), useRaster = TRUE, axes = FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
    legend(grconvertX(0.5, "device"), grconvertY(1, "device"),
           c(min(m), round(max(m), digits = 3)),
           fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
    par(old.par)
    return()
}

#### Threshold plot ################################################################
#' Threshold plot
#'
#' Plot cluster membership values and thresholding
#'
#' @param bce A BiclusterExperiment
#' @param bcs A BiclusterStrategy object, or the name or index of one of the
#'  BiclusterStrategy objects in \code{bce@@strategies}
#' @param type Specifies whether to show thresholding for features or samples
#' @param bicluster The character name of a bicluster in \code{bcs}. See
#'  \code{\link{bcNames}}.
#' @param ... Optional arguments \code{thresholds}, \code{ordering}, and
#'  \code{xlabs}, described below.
#' @param thresholds If not NULL, must be a numeric vector of theshold(s) to
#'  replace the algorithmically chosen thresholds contained in \code{bcs}.
#' @param ordering Specifies the ordering of features or samples along the
#'  x-axis. If "input", the input ordering is used. If "distance", ordering is
#'  based on Euclidean distance calculated from \code{bce@@assayData}. If
#'  "cluster", ordering is based on the bicluster membership of the samples or
#'  features.
#' @param xlabs Should sample or feature labels be displayed?
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' bcNames(bcs)
#' plotThreshold(bce, bcs, type = "sample", bicluster = "Bicluster.1")
#'
#' @return Invisibly returns NULL
#' @export
setGeneric("plotThreshold", signature = c("bce", "bcs"),
           function(bce, bcs, type, bicluster, ...) {
               standardGeneric("plotThreshold")
           })
#' @describeIn plotThreshold Plots thresholding performed in the
#'   \code{\link{BiclusterStrategy-class}} object named \code{bcs} in
#'   \code{bce@@strategies}.
setMethod("plotThreshold",
          signature(bce = "BiclusterExperiment", bcs = "character"),
          function(bce, bcs, type, bicluster, ...) {
              plotThreshold(bce = bce, bcs = getStrat(bce, bcs), type, bicluster, ...)
          }
)
#' @describeIn plotThreshold Plots thresholding performed in the
#'   \code{\link{BiclusterStrategy-class}} object at the index specified by
#'   \code{bcs}
setMethod("plotThreshold",
          signature(bce = "BiclusterExperiment", bcs = "numeric"),
          function(bce, bcs, type, bicluster, ...) {
              plotThreshold(bce, getStrat(bce, bcs), type, bicluster, ...)
          }
)
#' @describeIn plotThreshold Plots thresholding performed in \code{bcs}.
setMethod(
    "plotThreshold", signature(bce = "BiclusterExperiment", bcs = "BiclusterStrategy"),
    function(bce, bcs, type = c("feature", "sample"), bicluster = "Bicluster.1",
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
        data <- if(type == "feature") {
            logicalMatrix2Numeric(featureFactor(bcs, allBc = TRUE))
        } else {
            t(logicalMatrix2Numeric(sampleFactor(bcs, allBc = TRUE)))
        }
        data <- data[, bicluster, drop = TRUE]

        ordering <- match.arg(ordering)
        # Plot
        xs <- seq_along(data)
        xlab <- type
        ylab = "value"

        if(ordering == "input") {
            ord <- xs
        } else if (ordering == "distance") {
            m <- if(type == "sample") as.matrix(bce) else t(as.matrix(bce))
            ord <- hclust(dist(m))$order
        } else if (ordering == "cluster") {
            clusFunction <- if(type == "sample") {
                function(bcs) t(clusteredSamples(bcs))
            } else {
                function(bcs) clusteredFeatures(bcs)
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
        return()
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
        var <- ncol(as.matrix(bce))
        xs <- strat@factors@fit@H[1, ]
        var1 <- round(sd(xs) ^ 2 / var, 1)
        ys <- strat@factors@fit@H[2, ]
        var2 <- round(sd(ys) ^ 2 / var, 1)
    } else {
        m <- t(as.matrix(bce)) # pca on columns
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

# Create annotations dataframe for heatmaps
#
# Checks phenoLabels and biclustLabels for validity and then joins the
# corresponding data extracted from the provided BiclusterExperiment.
#
# Annotation tracks are converted to named dataframe columns.
#
# @param x a BiclusterExperiment
# @param names the names of entities annotated- required
# @param strategy index or name of a strategy in x. Required whenever
#   length(biclustLabels) > 0
# @param phenoLabels any phenotype labels in x
# @param biclustLabels any bicluster labels in x
createAnnots <-
    function(x, type = c("sample", "feature"),
             strategy = "",
             phenoLabels = c(),
             biclustLabels = c()) {
        annots <- NA
        # Process phenotype labels if type = samples
        if(type == "sample") {

            phenoLabels <- validatePhenoNames(phenoLabels, x)
            if (length(phenoLabels) > 0) {
                phData <-
                    as.data.frame(Biobase::pData(Biobase::phenoData(x))[, phenoLabels])
                colnames(phData) <- phenoLabels
                rownames(phData) <- rownames(Biobase::phenoData(x))
            }
        }

        # Process bicluster labels
        if (length(strategy) > 0 && length(biclustLabels) == 0) {
        } else if (length(biclustLabels) > 0) {
            if(!inherits(strategy, "BiclusterStrategy")) {
                if(nchar(strategy) == 0) {
                    stop("Argument 'strategy' missing", call. = FALSE)
                }
                strategy <- getStrat(x, strategy) # get BiclusterStrategy
            }
            biclustLabels <- validateBiclustNames(biclustLabels, strategy)
            predData <- as.data.frame(
                if(type == "sample") {
                    clusteredSamples(strategy)[biclustLabels, ]
                } else {
                    clusteredFeatures(strategy)[, biclustLabels]
                }
            )

            # bicluster annotations should be discrete
            predData[] <- as.data.frame(lapply(predData, as.factor))
            colnames(predData) <- biclustLabels
        }

        # Concatenate phenotype and bicluster labels if both requested
        annots <-
            if (length(phenoLabels) > 0 && length(biclustLabels) > 0) {
                cbind(phData, predData)
            } else if (length(phenoLabels) > 0) {
                phData
            }
        else if (length(biclustLabels) > 0) {
            predData
        }

        # pheatmap throws cryptic error without rownames
        if(inherits(annots, "data.frame")) {
            row.names(annots) <- if(type == "sample") {
                colnames(x)
            } else {
                rownames(x)
            }
        }

        return(annots)
    }

