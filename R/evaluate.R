# example usage
f <- function() {
  gds <- "GDS181" # 3/4 of genes have at least one NA
  dir.create(paste0("data/", gds))
  gds <- GEOquery::getGEO(gds, destdir = paste0("data/", gds))
  GEOquery::Meta(gds)$platform
  GEOquery::Table(gds)[1:10, 1:6]
  eSet <- GEOquery::GDS2eSet(gds)
  rm(gds)
  bce <- as(eSet, "BiclusterExperiment")
  bce <- clean(bce, 0.5) # Make sure to warn the user that their biclusterstrategies will be invalidated.
  bce <- addStrat(bce, bcs = BiclusterStrategy(m = t(as.matrix(bce)), k = 4, method = "nipals"))
}
# 


#' #' Clustering error
#' #' @param found a list of found biclusters, which consist of a tuple containing
#' #'   row range (as tuple) and column range (tuple)
#' scoreCe <- function(found, reference, mat) {
#'   ndot <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#'   lapply(found, function(bicluster) {
#'     ndot[bicluster[1][1]:bicluster[1][2]] <- 
#'       ndot[bicluster[1][1]:bicluster[1][2]] + 1
#'   })
#'   
#'   n <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
#'   lapply(found, function(bicluster) {
#'     n[bicluster[1][1]:bicluster[1][2]] <- 
#'       n[bicluster[1][1]:bicluster[1][2]] + 1
#'   })
#'   
#'   union <- sum(unlist(lapply(seq_len(len(mat)), function(i) {
#'     max(ndot[i], n[i])
#'   })))
#'   
#'   confusion <- matrix(0, length(reference), length(found))
#'   lapply(seq_len(length(reference)), function(i) {
#'     lapply(seq_len(length(found)), function(j) {
#'       confusion[i, j] <- union(found, reference)
#'     }
#'   }
#'   
#'   dmax <- sum(unlist(lapply(seq_len(min(length(found), length(reference))),
#'                             function(x) {}
#'   )))
#'                            
#' }
