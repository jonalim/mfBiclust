# example usage

evaluateGDS <- function() {
  # define list of files and known number biclusters
  gds <- c("GSE1", "GSE2223")
  gds <- "GSE3933"
  gds <- "GSE68907"
  #data/GSE1/GPL7.soft # grep "melanoma", others in source_name_ch1 # 19, 19
  # data/GSE2223/GPL1833.soft # grep case-sensitive BRAIN GBM O others in source_name_ch2 # 5, 31, 14, 5
  # data/GSE17025/GSE17025_series_matrix.txt.gz # grep EE, PS, and NL in sourece_name_ch1
  # data/GSE3726/GSE3726_series_matrix.txt.gz # grep B, C in title # 62, 42
  # "data/GSE4045/GSE4045_series_matrix.txt.gz" # grep serrated in description # 8, 29
  # data/GSE82009/GPL8300.soft # characteristics_ch1 can be read as factor # 	14,7,14,15
  # Ramaswamy multicancer # get labels 1:14 from labels files # tab-delimited
  # data/GSE68895/GPL80.soft # 2 biclusters. either grep follicular in source_name_ch1 or -- in characteristics_ch1.2
  
  # GSE3398 needs individual samples assembled (Garber adenocarcinomas)
  # "data/GSE3933/GSE3933-GPL3044_series_matrix.txt.gz" 43008 features (not processed into genes)

  # 
  oracle <- 
  # inside loop: fetch file.
    # perform biclustering
  #evaluate
  
}
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
biclusterGDS <- function(bce, k = 5) {
  biclusterTranscriptomicsHelper(bce, k, 0)
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
