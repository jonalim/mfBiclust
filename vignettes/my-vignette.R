## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options("max.print" = 80)

## ------------------------------------------------------------------------
library("mfBiclust")
set.seed(12345) # set seed for reproducibility
class(yeast_benchmark[[1]])
dim(yeast_benchmark[[1]])

dim(cancer_benchmark[[1]]$data)
cancer_benchmark[[1]]$labels

## ---- eval = FALSE-------------------------------------------------------
#  biclusterGUI()
#  biclusterGUI(yeast_benchmark[[1]])

## ------------------------------------------------------------------------
bce <- BiclusterExperiment(yeast_benchmark[[1]])
bce

## ------------------------------------------------------------------------
bce <- addStrat(bce, k = 3, method = "als-nmf")

## ------------------------------------------------------------------------
bcs <- getStrat(bce, 1)
# Equivalent to getStrat(bce, "ALS-NMF | Otsu | 3")

sampleFactor(bcs)
featureFactor(bcs)

## ---- fig.show = "hold"--------------------------------------------------
library("pheatmap")
# Which samples associate with which bicluster?
factorHeatmap(bce, bcs, type = "sample", colNames = TRUE, ordering = "cluster")

## ------------------------------------------------------------------------
# Which features are strongly present in which bicluster?
factorHeatmap(bce, bcs, type = "feature", colNames = FALSE, ordering = "cluster")

## ------------------------------------------------------------------------
bcNames(bcs)
plotThreshold(bce = bce, bcs = bcs, type = "sample", bicluster = "Bicluster.1",
                  ordering = "cluster", xlabs = TRUE)
# ordering = "input" plots samples in the same order as their respective columns
# in the BiclusterExperiment

## ------------------------------------------------------------------------
plotThreshold(bce = bce, bcs = bcs, type = "feature", bicluster = "Bicluster.1",
                  ordering = "cluster", xlabs = FALSE)

## ------------------------------------------------------------------------
clusteredSamples(bcs)
clusteredFeatures(bcs)

## ------------------------------------------------------------------------
# Samples and features in the first bicluster
names(which(clusteredSamples(bcs)[1, ]))
names(which(clusteredFeatures(bcs)[, 1]))

# Samples and features in the second bicluster
 names(which(clusteredSamples(bcs)[2, ]))
names(which(clusteredFeatures(bcs)[, 2]))

# and so on

## ------------------------------------------------------------------------
auto_bcv(Y = as.matrix(bce), bestOnly = FALSE)
# To omit the counts table, omit `bestOnly` argument

## ------------------------------------------------------------------------
bce <- addStrat(bce, k = 3, method = "snmf", beta = 3, eta = 3)
factorHeatmap(bce, getStrat(bce, 2), type = "sample", colNames = TRUE, ordering = "cluster")

## ------------------------------------------------------------------------
factorHeatmap(bce, getStrat(bce, 2), type = "feature", colNames = FALSE, ordering = "cluster")

## ------------------------------------------------------------------------
bce <- addStrat(bce, k = 3, method = "svd-pca")
factorHeatmap(bce, getStrat(bce, 3), type = "sample", colNames = TRUE, ordering = "cluster")

## ------------------------------------------------------------------------
factorHeatmap(bce, getStrat(bce, 3), type = "feature", colNames = FALSE, ordering = "cluster")

