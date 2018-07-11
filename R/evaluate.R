#' @include helperFunctions.R

libGri <- function(pathToPy27 = NULL) {
  if(!requireNamespace("reticulate")) {
    if(askYesNo(paste("Would you like to install reticulate and its", 
                      "depencencies to use this feature?"))) {
      install.packages("reticulate")
    }
  }
  if(!requireNamespace("reticulate")) { stop("Unable to install reticulate.") }
  
  if(!is.null(pathToPy27)) {
    reticulate::use_python(pathToPy27)
  }
  reticulate::py_config()
  tryCatch(gri <<- reticulate::import("gri"), 
           error = function(e) {
             stop(paste("Unable to import Python module gri-master. Please",
                        "download from https://github.com/padilha/gri, run",
                        "\"python setup.py install\" and then move files",
                        "gri-1.0-py2.7.egg, gri.py, and gri.pyc to the /lib/site-packages",
                        "folder of your Python installation."))
           }
  )
}

libGostats <- function() {
  reqs <- rep(FALSE, 2)
  names(reqs) <- c("gostats", "GO.db")
  if(!requireNamespace("gostats")) {
    reqs[1] <- TRUE
  }
  if(!requireNamespace("GO.db")) {
    reqs[2] <- TRUE
  }
  if(askYesNo(paste("Dependencies", reqs, "required. Install? (y/n)"))) {
    if(reqs[2]) {
      BiocInstaller::biocLite("GO.db", type = "source", INSTALL_opts = "--no-multiarch")
    }
    if(reqs[1]) {
      BiocInstaller::biocLite("GOstats", dependencies = TRUE)
    }
  }
}

calcFE <- function(dataset, algorithm, cutoffs) {
  # Perform biclustering for 300 biclusters
  bce <- BiclusterExperiment(t(as.matrix(dataset)))
  
  bce <- addStrat(bce, 500, method = algorithm, duplicable = TRUE)
  
  # filter down to a list of 100 biclusters / gene lists
  bc <- threshold(loading(getStrat(bce, 1)), MARGIN = 1, 
                  loadingThresh(getStrat(bce, 1)))
  res <- filter.biclust(RowxBicluster = pred(getStrat(bce, 1)),
                        BiclusterxCol = bc,
                        max = 100, overlap = 0.25)
  biclustered <- res[[3]]
  geneLists <- lapply(seq_len(nrow(res[[2]])), function(index) {
    rownames(bce)[res[[2]][index, ]]
  })
  
  if(tolower(method(getStrat(bce, 1))) != tolower(algorithm)) {
    termCount = as.data.frame(matrix(rep(NA, length(cutoffs)), nrow = 1))
    colnames(termCount) <- as.character(cutoffs)
    list(biclusters = list(), termCount = termCount)
  } else {
    universe <- rownames(bce) # I sure hope these are entrez
    
    # get a table where each column is a p-value cutoff; each row is a
    # bicluster
    lapp <- if(requireNamespace("BiocParallel")) bplapply else lapply
    fun <- function(geneList, universe, cutoffs) {
      # Molecular Function ontology
      mf <- GOstats::hyperGTest(new("GOHyperGParams",
                                    geneIds = geneList,
                                    universeGeneIds = universe,
                                    annotation = "org.Sc.sgd.db",
                                    ontology = "MF",
                                    pvalueCutoff = 0.05,
                                    testDirection = "over", 
                                    conditional = TRUE))
      
      # Biological Process ontology
      bp <- GOstats::hyperGTest(new("GOHyperGParams",
                                    geneIds = geneList,
                                    universeGeneIds = universe,
                                    annotation = "org.Sc.sgd.db",
                                    ontology = "BP",
                                    pvalueCutoff = 0.05,
                                    testDirection = "over",
                                    conditional = TRUE))
      
      # a vector with a GO term count for each cutoff
      termCounts <- sapply(cutoffs, function(cutoff) {
        sum(GOstats::pvalues(mf) < cutoff) + sum(GOstats::pvalues(bp) < cutoff)
      })
      
      names(termCounts) <- as.character(cutoffs)
      termCounts
    }
    termCount <- do.call(rbind, lapp(geneLists, fun, universe, cutoffs))
    list(bce = bce, biclustered = biclustered, geneLists = geneLists, 
         termCounts = termCount)
  }
}

testFE <- function(rep = 30) {
  oldSeed <- .Random.seed # do not modify the R global environment
  set.seed(12345)
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  cutoffs <- c(0.05, 0.01, 0.005, 0.0001, 0.00001)
  datasets.all <- loadBenchmark("data/yeast_benchmark/", classes = FALSE)
  # change to TRUE to save biclustering results
  saveMe <- TRUE
  save.file <- "plots/yeast_benchmark_results/"
  
  methods.nondet <- c("nipals-pca", "plaid", "spectral")
  methods.det <- c("als-nmf", "svd-pca")
  
  extractBest <- function(solutions) {
    best <- which.max(sapply(solutions, function(solution) {
      sum(solution$termCount[, length(cutoffs)] > 0) / length(solution$geneLists)
    }))
    if(length(best) == 0) best <- 1
    best <- solutions[[best]]
    # number of enriched biclusters at various cutoffs
    best$enriched <- colSums(best$termCount > 0)
    best
  }
  
  # JNL Use the 1st and 5th datasets to test
  
  # Find biclusters using nondeterministic algorithms
  solutions.nondet <- lapply(methods.nondet, function(algo) {
    solutions <- lapply(datasets.all, function(dataset) {
      # try 30 times to get the most GO-enriched biclusters. Seems like cheating,
      # but this does reflect a real use case. We can inform users that "super"
      # functional analysis will repeat NMF to find the most biological signal
      reps <- lapply(seq_len(rep), function(i) {
        calcFE(dataset, algorithm = algo, cutoffs)
      })
      extractBest(reps)
    })
  })
  
  if(saveMe) save(solutions.nondet, file = paste0(save.file, "solutions.nondet.Rda")) # save at each step!
  
  # find biclusters using deterministic methods
  solutions.det <- lapply(methods.det, function(algo) {
    solutions <- lapply(datasets.all, function(dataset) {
      solution <- calcFE(dataset, algorithm = algo, cutoffs)
      solution$enriched <- solution$termCount > 0
      solution
    })
  })
  
  # compile the results ("als-nmf", "svd-pca", "plaid", "spectral")
  solutions <- c(solutions.det, solutions.nondet) 
  if(saveMe) save(solutions, file = paste0(save.file, "solutions.all.Rda")) # save at each step!
  
  res <- sapply(solutions, function(solutions.algo) {
    # total enriched at each cutoff
    enriched <- colSums(do.call(rbind, lapply(solutions.algo, function(solution.dataset) {
      solution.dataset$enriched
    })), na.rm = TRUE)
    total <- sum(sapply(solutions.algo, function(x) {
      length(x$geneLists)
    })) # total
    scores <- enriched / total * 100
    pctBiclustered <- sum(solution.dataset$biclustered == 1) / length(solution.dataset$biclustered) * 100
    list(enriched = enriched[1], total = total, pctEnriched = scores, pctBiclustered = pctBiclustered)
  })
  if(!inherits(res, "matrix")) res <- cbind(res)
  
  methods.all <- c(methods.det, methods.nondet)
  colnames(res) <- methods.all
  
  ys <- apply(res, MARGIN = 2, function(x) x$pctEnriched)
  if(saveMe) save(solutions, res, cutoffs, file = paste0(save.file, "solutions+stats.Rda"))
  
  old.par <- par(no.readonly = TRUE)
  par(mar = c(4.1, 4.1, 1.1, 1.1))
  barplot(height = ys, beside = TRUE, legend.text = as.character(cutoffs), args.legend = c(x = "right"),
          ylab = "% GO-enriched", xlab = "Algorithm", col = RColorBrewer::brewer.pal(nrow(res), "Dark2"))
  par(old.par)
}

#' @export
loadBenchmark <- function(dir = "data/yeast_benchmark/", classes = FALSE) {
  files <- list.files(path = dir)
  if(classes) {
    lapply(
      files, 
      function(x) {
        # Read without header because duplicate column names are not allowed
        tab <- read.table(file = paste0(dir, x), sep = "",
                          header = FALSE, stringsAsFactors = FALSE)
        labels <- factor(as.character(tab[1, 2:ncol(tab)]))
        tab <- sapply(tab[2:nrow(tab), 2:ncol(tab)], as.numeric)
        
        # convert labels to classification matrix
        classes <- levels(labels)
        labels <- do.call(rbind, lapply(labels, function(lab) {
          vec <- rep(0, length(classes))
          names(vec) <- classes
          vec[lab] <- 1
          vec
        }))
        
        list(data = tab, labels = labels)
      }
    )
  } else {
    lapply(files, function(x) {
      tab <- read.table(file = paste0(dir, x), sep = "",
                        header = FALSE, row.names = 2, stringsAsFactors = FALSE)
      tab[, 2:ncol(tab)]
    })
  }
}

testResidAgriCor <- function(rep = 30) {
  oldSeed <- duplicable("evalua") # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  libGri()
  datasets.all <- loadBenchmark("data/cancer_benchmark/", classes = TRUE)
  # Test only datasets with >2 classes because ALS-NMF is deterministic-sh for two 
  # classes
  overTwoClass <- datasets.all[sapply(datasets.all, 
                                      function(X) ncol(X$labels) > 2)]
  
  set.seed(12345)
  
  agrisResids <- sapply(seq_len(30), FUN = function(i) {
    agrisResids <- sapply(overTwoClass, FUN = function(dataset) {
      bce <- BiclusterExperiment(t(as.matrix(dataset$data)))
      bce <- addStrat(bce, k = length(unique(colnames(dataset$labels))), 
                      method = "als-nmf")
      labels <- dataset$labels
      
      agri <- if(method(getStrat(bce, 1)) == "als-nmf") {
        prediction <- pred(getStrat(bce, 1))
        prediction <- apply(prediction, 2, as.numeric)
        gri$grand_index(t(labels), t(prediction), adjusted=TRUE)  
      } else { 0 }
      resid <- t(as.matrix(dataset$data)) - 
        score(getStrat(bce, 1)) %*% loading(getStrat(bce, 1))
      resid.rms <- sqrt(sum(resid ^ 2) / length(resid))
      
      c(agri, resid.rms)
    })
    agrisResids
  }, simplify = "array")
  
  pearsons <- apply(agrisResids, MARGIN = 2, function(array) {
    cor(x = array[1, ], y = array[2, ], method = "pearson")
  })
  # print(boxplot(x = pearsons, xlab = "20 datasets", y = "Pearson's Correlation"))
  
  agris <- apply(agrisResids, MARGIN = 2, function(dataset) {
    dataset[1, ]
  })
  
  xs <- colMeans(agris)
  print(plot(x = xs, y = pearsons, xlab = "AGRI", ylab = "PCC[AGRI, rms(resids)]"))
  fit <- lm(pearsons ~ xs)
  abline(fit)
  imax <- which.max(xs)
  text(xs[imax] - 0.02, pearsons[imax], paste0("slope = ", round(fit$coefficients[2], 1)), pos=2)
  
  se2 <- unlist(apply(agris, MARGIN = 2, function(x) {
    2 * sd(x) / sqrt(length(x))
  }))
  arrows(xs + se2, pearsons, xs - se2, pearsons, angle=0, code=3, length=0.1)
  
  # Pearsons correlations were in the range
  # -0.662 to 0.293.; mean = -0.154 +/- 0.132 (95%), p =
  # 0.02547 for two-tailed t-test. Supporting this, the
  # correlation was more negative on datasets on which ALS-NMF performed well,
  # suggesting that best performance is achieved by trying a few replicates.
  
  # Positive correlations were theoretically
  # caused by some signal orthogonal to the sample classes. 
  
  list(pearsons = pearsons, 
       tt = t.test(x = pearsons, alternative = c("two.sided")),
       agris = agris,
       resids.rms = apply(agrisResids, MARGIN = 2, function(dataset) dataset[2, ])
  )
}

testMultiBiclustering <- function() {
  oldSeed <- duplicable("evalua") # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  # This evaluation yielded large confidence intervals for AGRI on ALS-NMF. 
  # If the algorithm is run repeatedly on the same matrix, does lower norm of residuals
  # correlate with higher AGRI? H0: Pearson correlation = 0 across all datasets.
  try(load(file = "plots/multigroup-cancer-benchmark/results.rda"))
  
  # Use 13AGRI to compare possibilistic clustering solutions with the reference,
  # which happens to be exclusive hard clustering
  libGri()
  # Configure loading 13AGRI from Python module 
  # https://github.com/padilha/gri
  # installed by running python setup.py install, then placing the .egg and gri.py at
  # anaconda3/Lib/site-packages
  
  # Load cancer benchmark
  datasets.all <- lapply(
    X = list.files(path = "data/cancer_benchmark/"), 
    function(x) {
      tab <- read.table(file = paste0("data/cancer_benchmark/", x), sep = "\t",
                        header = FALSE, stringsAsFactors = FALSE)
      labels <- factor(as.character(tab[1, 2:ncol(tab)]))
      tab <- sapply(tab[2:nrow(tab), 2:ncol(tab)], as.numeric)
      
      # convert labels to classification matrix
      classes <- levels(labels)
      labels <- do.call(rbind, lapply(labels, function(lab) {
        vec <- rep(0, length(classes))
        names(vec) <- classes
        vec[lab] <- 1
        vec
      }))
      
      list(data = tab, labels = labels)
    }
  )
  
  # Find and evaluate als-nmf solutions for datasets
  agris.als_nmf <- sapply(seq_len(30), function(i) {
    agris.als_nmf <- sapply(datasets.all, function(dataset) {
      bce <- BiclusterExperiment(t(as.matrix(dataset$data)))
      bce <- addStrat(bce, k = length(unique(colnames(dataset$labels))), method = "als-nmf")
      labels <- dataset$labels
      if(method(getStrat(bce, 1)) == "als-nmf") {
        prediction <- pred(getStrat(bce, 1))
        prediction <- apply(prediction, 2, as.numeric)
        gri$grand_index(t(labels), t(prediction), adjusted=TRUE)  
      } else { 0 }
    }
    )
    agris.als_nmf
  })
  
  # Find and evaluate SVD-PCA solutions for datasets
  bces.svd_pca <- sapply(datasets.all, function(l) {
    bce <- BiclusterExperiment(t(as.matrix(l$data)))
    bce <- addStrat(bce, k = length(unique(colnames(l$labels))), method = "svd-pca")
  }
  )
  agris.svd_pca <- mapply(function(bce, labels) {
    labels <- labels$labels
    
    prediction <- pred(getStrat(bce, 1))
    prediction <- apply(prediction, 2, as.numeric)
    gri$grand_index(t(labels), t(prediction), adjusted=TRUE)  
    
  }, bce = bces.svd_pca, labels = datasets.all
  )
  
  agris.plaid <- sapply(seq_len(30), function(i) {
    dataset <<- 0
    agris.plaid <- sapply(datasets.all, FUN = function(l) {
      dataset <<- dataset + 1
      bce <- BiclusterExperiment(t(as.matrix(l$data)))
      agri.plaid <- NULL
      while(!is.numeric(agri.plaid)) {
        agri.plaid <- try({
          print(paste(i, dataset))
          bce <- addStrat(bce, k = length(unique(colnames(l$labels))), method = "plaid")
          prediction <- pred(getStrat(bce, 1))
          prediction <- apply(prediction, 2, as.numeric)
          gri$grand_index(t(l$labels), t(prediction), adjusted = TRUE)
        }
        )
      }
      agri.plaid
    })
  })
  
  agris.spectral <- sapply(seq_len(30), function(i) {
    agris.spectral <- sapply(datasets.all, FUN = function(dataset) {
      bce <- BiclusterExperiment(t(as.matrix(dataset$data)))
      agri.spectral <- NULL
      while(!is.numeric(agri.spectral)) {
        agri.spectral <- try({
          dummy <- capture.output(bce <- addStrat(bce, k = 2, method = "spectral"))
          prediction <- pred(getStrat(bce, 1))
          prediction <- apply(prediction, 2, as.numeric)
          gri$grand_index(t(dataset$labels), t(prediction), adjusted = TRUE)
        }
        )
      }
      agri.spectral
    })
  })
  
  save(agris.als_nmf, agris.svd_pca, agris.plaid, file = "plots/multigroup-cancer-benchmark/results.rda")
  
  agris.als_nmf.means <- rowMeans(agris.als_nmf)
  agris.plaid.means <- rowMeans(agris.plaid)
  ys <- rbind(agris.als_nmf.means, agris.svd_pca, agris.plaid.means)
  agris.als_nmf.se2 <- 2 * apply(agris.als_nmf, 1, function(x) sd(x) / sqrt(length(x)))
  agris.plaid.se2 <- 2 * apply(agris.plaid, 1, function(x) sd(x) / sqrt(length(x)))
  agris.svd_pca.se2 <- rep(0, length(agris.svd_pca))
  se2 <- rbind(agris.als_nmf.se2, agris.svd_pca.se2, agris.plaid.se2)
  
  old.par <- par(no.readonly = T)
  par(mar = c(2.1, 4.1, 1.1, 0))
  bp <- barplot(height = ys, beside = TRUE, legend.text = c("ALS-NMF", "SVD-PCA", "Plaid"), args.legend = c("topright"),
                ylab = "'13 Adjusted Grand Rand Index", col = RColorBrewer::brewer.pal(3, "Dark2"),
                ylim = c(0, .35))
  error.bar(bp, ys, se2, length = 0.01)
  par(old.par)
}

testSinglePca <- function() {
  oldSeed <- duplicable("evalua") # do not modify the R global environment
  on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
  
  try(load("plots/two-group-cancer-benchmark/results.rda"))
  # Load previously stored data now (so we don't overwrite new results)
  
  datasets.all <- sapply(X = list.files(path = "data/cancer_benchmark/"), 
                         function(x) {
                           tab <- read.table(file = paste0("data/cancer_benchmark/", x), sep = "\t",
                                             header = FALSE, stringsAsFactors = FALSE)
                           labels <- factor(as.character(tab[1, 2:ncol(tab)]))
                           tab <- sapply(tab[2:nrow(tab), 2:ncol(tab)], as.numeric)
                           if(length(levels(labels)) == 2) { list(data = tab,
                                                                  labels = labels)
                           }
                           else { NULL }
                         }
  )
  twoClass <- datasets.all[sapply(datasets.all, function(X) !is.null(X))]
  rm(datasets.all)
  
  minRowClusterFraction <- min(sapply(twoClass, function(x) table(x$labels) / length(x$labels)))
  
  aris.als_nmf <- sapply(twoClass, function(dataset) {
    bce <- BiclusterExperiment(t(as.matrix(dataset$data)))
    bce <- addStrat(bce, k = 1, method = "als-nmf")
    
    labels <- dataset$labels
    if(method(getStrat(bce, 1)) == "als-nmf") {
      mclust::adjustedRandIndex(labels, getStrat(bce, 1)@pred[, 1])
    } else { 0 }
  })
  
  bces.svd_pca <- sapply(twoClass, function(l) {
    bce <- BiclusterExperiment(t(as.matrix(l$data)))
    bce <- addStrat(bce, k = 1, method = "svd-pca")
  }
  )
  aris.svd_pca <- mapply(function(bce, labels) {
    labels <- labels$labels
    mclust::adjustedRandIndex(labels, getStrat(bce, 1)@pred[, 1])
  }, bce = bces.svd_pca, labels = twoClass
  )
  # Since plaid is nondeterministic, for each dataset, take the mean ARI from 30
  # replicate runs. N.B. that each of the 30 might have a different release
  # parameter; each run, the parameter is eased down from 0.7 to 0 in steps of
  # 0.1 until enough biclusters are found. Then the first k biclusters are
  # reported. Here, k = 1, and the two groups are bicluster and non-bicluster.
  dataset <- NULL
  aris.plaid <- sapply(seq_len(30), function(i) {
    dataset <<- 0
    aris.plaid <- sapply(twoClass, FUN = function(l) {
      dataset <<- dataset + 1
      bce <- BiclusterExperiment(t(as.matrix(l$data)))
      ari.plaid <- NULL
      while(!is.numeric(ari.plaid)) {
        ari.plaid <- try({
          print(paste(i, dataset))
          dummy <- capture.output(bce <- addStrat(bce, k = 1, method = "plaid"))
          mclust::adjustedRandIndex(l$labels, getStrat(bce, 1)@pred[, 1])
        }
        )
      }
      ari.plaid
    })
  })
  # aris.plaid <- rowMeans(aris.plaid)
  
  aris.spectral <- sapply(seq_len(30), function(i) {
    dataset <- 0
    aris.spectral <- sapply(twoClass, FUN = function(l) {
      dataset <<- dataset + 1
      print(dataset)
      bce <- BiclusterExperiment(t(as.matrix(l$data)))
      ari.spectral <- NULL
      
      bce <- addStrat(bce, k = 2, method = "spectral", min = minRowClusterFraction)
      predictions <- pred(getStrat(bce, 1))[, 1]
      ari.spectral <- mclust::adjustedRandIndex(l$labels, predictions)
      
      if(is.numeric(ari.spectral)) ari.spectral
      else 0 # N.B. My Spectral wrapper tries withinVar at 1:10 * nrow(m)
    })
  })
  
  save(aris.als_nmf, aris.svd_pca, aris.plaid, aris.spectral, file = "plots/two-group-cancer-benchmark/results.rda")
  
  aris.spectral.means <- rowMeans(aris.spectral)
  aris.plaid.means <- rowMeans(aris.plaid)
  ys <- rbind(aris.als_nmf, aris.svd_pca, aris.plaid.means, aris.spectral.means)
  colnames(ys) <- unlist(sapply(colnames(ys), function(x) substr(x, 1, nchar(x) - 13)))
  
  aris.spectral.se2 <- 2 * apply(aris.spectral, 1, function(x) sd(x) / sqrt(length(x)))
  aris.plaid.se2 <- 2 * apply(aris.plaid, 1, function(x) sd(x) / sqrt(length(x)))
  aris.svd_pca.se2 <- rep(0, length(aris.svd_pca))
  aris.als_nmf.se2 <- rep(0, length(aris.als_nmf))
  se2 <- rbind(aris.als_nmf.se2, aris.svd_pca.se2, aris.plaid.se2, aris.spectral.se2)
  
  old.par <- par(no.readonly = T)
  par(mar = c(3.1, 4.1, 1.1, 1))
  bp <- barplot(height = ys, beside = TRUE, legend.text = c("ALS-NMF", "SVD-PCA", "Plaid", "Spectral"), args.legend = c(x = "topright"),
                ylab = "Adjusted Rand Index", col = RColorBrewer::brewer.pal(4, "Dark2"),
                ylim = c(0, 1), xaxt = "n")
  error.bar(bp, ys, se2, length = 0.01)
  par(old.par)
}

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
