# example usage
geo2Bce <- function(gds = "gds181") {
  # 3/4 of genes have at least one NA
  dir.create(paste0("data/", gds))
  gds <- GEOquery::getGEO(gds, destdir = paste0("data/", gds))
  GEOquery::Meta(gds)$platform
  GEOquery::Table(gds)[1:10, 1:6]
  eSet <- GEOquery::GDS2eSet(gds)
  bce <- as(eSet, "BiclusterExperiment")
}

#' @export
biclusterTranscriptomics <- function(bce, maxK = 500) {
  biclusterTranscriptomicsHelper(bce, maxK, 0)
}

biclusterTranscriptomicsHelper <- function(bce, maxK, cleanParam = 0) {
  if(cleanParam > 0) { bce <- clean(bce, cleanParam) }
  
  tryCatch({
    addStrat(bce, bcs = BiclusterStrategy(m = t(as.matrix(bce)), k = maxK, method = "nipals"))
  }, error = function(e) {
    cleanParam <- cleanParam + (1 - cleanParam) / 2
    message(paste("Too many NA in the data. Cleaning with maxNAs at", 
                  cleanParam))
    biclusterTranscriptomicsHelper(bce, maxK, cleanParam)
  })
}

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
