#' @include helperFunctions.R
# example usage

#' @importFrom Biobase pData
testMicroarrays <- function() {
  # define list of files and known number biclusters
  
  #data/GSE1/GPL7.soft # grep "melanoma", others in source_name_ch1 # 19, 19
  gse1 <- GEOquery::getGEO(file = "data/annotated_microarrays/GSE1/GSE1_series_matrix.txt.gz")
  whichClustered <- read.csv(file = "data/annotated_microarrays/GSE1/which_clustered.csv", stringsAsFactors = FALSE)
  gse1.labels <- rep("unclustered", nrow(pData(gse1)))
  gse1.labels[which(sapply(pData(gse1)$title, function(x) {
    s <- substring(text = x, first = 7, last = 30)
    stringi::stri_replace_first(s, ".", regex = "-")
    }) %in% whichClustered[, 1])] <- "clustered"
  
  # data/GSE2223/GPL1833.soft # grep case-sensitive BRAIN GBM O others in source_name_ch2 # 4, 31, 14, 6
  gse2223 <- GEOquery::getGEO(file = "data/annotated_microarrays/GSE2223/GSE2223_series_matrix.txt.gz")
  gse2223.labels <- rep("Other", nrow(pData(gse2223)))
  gse2223.labels[grep(pattern = "BRAIN", x = pData(gse2223)$source_name_ch2, 
                     ignore.case = FALSE)] <- "BRAIN"
  gse2223.labels[grep(pattern = "GBM", x = pData(gse2223)$source_name_ch2, 
                     ignore.case = TRUE)] <- "GBM"
  gse2223.labels[grep(pattern = "GNN", x = pData(gse2223)$source_name_ch2, 
                      ignore.case = TRUE)] <- "GBM"
  gse2223.labels[grep(pattern = "O", x = pData(gse2223)$source_name_ch2, 
                      ignore.case = FALSE)] <- "O"
  
  # data/GSE17025/GSE17025_series_matrix.txt.gz # grep EE, PS, and NL in sourece_name_ch1
  # 79, 12, 12
  gse17025 <- GEOquery::getGEO(file = "data/annotated_microarrays/GSE17025/GSE17025_series_matrix.txt.gz")
  gse17025.labels <- rep(0, nrow(pData(gse17025)))
  gse17025.labels[grep(pattern = "EE", x = pData(gse17025)$source_name_ch1, 
                      ignore.case = FALSE)] <- "EE"
  gse17025.labels[grep(pattern = "PS", x = pData(gse17025)$source_name_ch1, 
                      ignore.case = TRUE)] <- "PS"
  gse17025.labels[grep(pattern = "NL", x = pData(gse17025)$source_name_ch1, 
                      ignore.case = TRUE)] <- "NL"
  
  exprSets <- list(gse1, gse17025, gse2223)
  labels <- list(gse1.labels, gse17025.labels, gse2223.labels)
  # do this for all of the 2-group microarrays...maybe even the 4-group ones?
  mapply(function(es, labels) {
    bce <- as(es, "BiclusterExperiment")
    k.oracle <- length(labels)
    addStrat(bce, BiclusterStrategy(t(as.matrix(bce)), k.oracle, method = "nipals"))
    adjustedRandIndex(gse1.labels, getStrat(bce, 1)@pred[, 1])
    },
    es = exprSets, labels = labels)
  # data1/GSE3726/GSE3726_series_matrix.txt.gz # grep B, C in title # 62, 42
  # "data/GSE4045/GSE4045_series_matrix.txt.gz" # grep serrated in description # 8, 29
  # data/GSE82009/GPL8300.soft # characteristics_ch1 can be read as factor # 	14,7,14,15
  # Ramaswamy multicancer # get labels 1:14 from labels files # tab-delimited
  # data/GSE68895/GPL80.soft # 2 biclusters. either grep follicular in source_name_ch1 or -- in characteristics_ch1.2
  
  # GSE3398 needs individual samples assembled (Garber adenocarcinomas)
  # "data/GSE3933/GSE3933-GPL3044_series_matrix.txt.gz" 43008 features (not processed into genes)

  # 
  # inside loop: fetch file.
    # perform biclustering
  #evaluate
  
}

#' Perform biclustering on microarray data
#' 
#' Returns a BiclusterExperiment object with Bicluster.X fields in phenoData.
#' Optionally, also a list containing bicluster assignments
#' @export
biclusterGDS <- function(bce, k = 5) {
  autoNipals() # FIXME
}

geo2Bce <- function(gds = "gds181") {
  # 3/4 of genes have at least one NA
  dir.create(paste0("data/", gds))
  gds <- GEOquery::getGEO(gds, destdir = paste0("data/", gds))
  GEOquery::Meta(gds)$platform
  GEOquery::Table(gse1)[1:10, 1:6]
  eSet <- GEOquery::GDS2eSet(gse1)
  bce <- as(gse1, "BiclusterExperiment")
}

setGeneric("performance", function(x) standardGeneric("performance"))

#' Automatically create a BiclusterStrategy for an existing 
#'
#' speed up by setting cleanParam to 0.5
#' 
#' @export
auto_bcs <- function(bce, k, cleanParam = 0) {
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
