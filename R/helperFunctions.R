biclusterNumber2scoreLoading <- function(biclusterNumber, m, k) {
  scores <- do.call(cbind, lapply(seq_len(k), function(i) {
    bicluster <- biclusterNumber[[i]]
    s <- rep(0, nrow(m))
    s[bicluster$Rows] <- 1
    s
  }))
  
  loadings <- do.call(rbind, lapply(seq_len(k), function(i) {
    bicluster <- biclusterNumber[[i]]
    l <- rep(0, ncol(m))
    l[bicluster$Cols] <- 1
    l
  }))
  
  list(scores, loadings)
}

#### capitalize ####
capitalize <- Vectorize(function(s) {
  # Use this function whenever names will be displayed in plots or gui
  if (s == "als-nmf") { "ALS-NMF" }
  if (s == "nipals-pca") { "NIPALS-PCA" }
  if (s == "svd-pca") { "SVD-PCA" }
  else { switch(s, snmf = "SNMF", pca = "PCA", otsu = "Otsu", 
                paste0(toupper(substring(s, 1,1)), substring(s, 2))) }
})

#' Create annotations dataframe for heatmaps
#'
#' Checks phenoLabels and biclustLabels for validity and then joins the
#' corresponding data extracted from the provided BiclusterExperiment.
#'
#' Annotation tracks are converted to named dataframe columns.
#'
#' @param x a BiclusterExperiment
#' @param names the names of entities annotated- required
#' @param strategy a strategy in x. Required whenever length(biclustLabels) > 0
#' @param phenoLabels any phenotype labels in x
#' @param biclustLabels any bicluster labels in x
createAnnots <-
  function(x,
           names,
           strategy = "",
           phenoLabels = c(),
           biclustLabels = c()) {
    annots <- NA
    # Process phenotype labels
    phenoLabels <- validatePhenoNames(phenoLabels, x)
    if (length(phenoLabels) > 0) {
      phData <-
        as.data.frame(Biobase::pData(Biobase::phenoData(x))[, phenoLabels])
      colnames(phData) <- phenoLabels
    }
    
    # Process bicluster labels
    if (length(strategy) > 0 && length(biclustLabels) == 0) {
      warning(
        paste0(
          "Since bicluster annotations were not ",
          "requested, the \"strategy\" argument will be ignored."
        )
      )
    } else if (length(biclustLabels) > 0) {
      validateStratName(strategy, x) # if no strategy, stop
      bcs <- getStrat(x, strategy) # get BiclusterStrategy
      biclustLabels <- validateBiclustNames(biclustLabels, bcs)
      predData <- as.data.frame(pred(bcs)[, biclustLabels])
      
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
    if (inherits(annots, "data.frame")) {
      row.names(annots) <- names
    }
    annots
  }

duplicable <- function(str) {
  if (!exists(".Random.seed", mode="numeric")) sample(NA)
  oldSeed <- .Random.seed
  newSeed <- strtoi(str, 35L)
  if(inherits(newSeed, "numeric")) set.seed(newSeed) else set.seed(12345)
  oldSeed
}

#' Error bars
#' 
#' Function adapted from James Holland Jones, 2009.
#' http://monkeysuncle.stanford.edu/?p=485
#' 
#' @param x the result of a call to barplot()
#' @param y the ys provided to barplot()
#' @param upper vector of lengths of upper error bars
#' @param lower vector of lengths of lower error bars; by default, same as 
#'   upper
#' @export
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  suppressWarnings(arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...))
}

filter.biclust <- function(RowxBicluster, BiclusterxCol, max = NULL, 
                           overlap = 0.25) {
  if(ncol(RowxBicluster) != nrow(BiclusterxCol)) { # validate
    stop(paste0("RowxBicluster must have the same number of columns as",
                "BiclusterxCol"))
  }
  
  k <- ncol(RowxBicluster)
  if(k == 1 || k == 0) {
    # no filtering needed
    list(RowxBicluster = RowxBicluster, 
         BiclusterxCol = BiclusterxCol)
  } else {
    
    # Create lists of rows and columns contained in biclusters
    BiclusterRows <- apply(RowxBicluster, MARGIN = 2, which)
    BiclusterCols <- apply(BiclusterxCol, MARGIN = 1, which)
    
    chosen <- rep(FALSE, each = k)
    pool <- rep(TRUE, each = k)
    sizes <- sapply(seq_len(k), function(biclus) {
      length(BiclusterRows[[biclus]]) * length(BiclusterCols[[biclus]])
    })
    names(sizes) <- seq_len(k)
    
    # exclude empty biclusters and whole-dataset biclusters
    pool[sizes == 0 | sizes == nrow(RowxBicluster) * ncol(BiclusterxCol)] <- FALSE
    
    # For each bicluster, calculate its overlap with all other biclusters
    overlaps <- lapply(seq_len(k), function(biclus1) {
      sapply(seq_len(k), function(biclus2) {
        intersection1 <- base::intersect(BiclusterRows[biclus1], BiclusterRows[biclus2])
        intersection1 <- if(length(intersection1) > 0) length(intersection1[[1]])
        else 0
        intersection2 <- base::intersect(BiclusterCols[biclus1], BiclusterCols[biclus2])
        intersection2 <- if(length(intersection2) > 0) length(intersection2[[1]])
        else 0
        intersection <- intersection1 * intersection2
        overlap <- intersection / sizes[biclus1]
      })
    })
    
    while(all(sum(chosen) < max) && sum(pool) > 0) {
      chooseMe <- as.numeric(names(which.max(sizes[which(pool)])))
      chosen[chooseMe] <- TRUE
      # Remove from the pool any biclusters heavily overlapping with chooseMe
      pool[overlaps[[chooseMe]] > overlap] <- FALSE
      pool[chooseMe] <- FALSE
    }
    list(RowxBicluster = RowxBicluster[, chosen], 
         BiclusterxCol = BiclusterxCol[chosen, ])
  }
}

is.wholenumber <-
  function(x, tol = sqrt(.Machine$double.eps)) {
    abs(x - round(x)) < tol
  }
# 
# autoNipals <- function(m, k, cleanParam = 0) {
#   cleanRes <- clean(m, cleanParam, index = TRUE)
#   mClean <- cleanRes$obj
#   indexRem <- cleanRes$indexRemaining
#   
#   tryCatch({
#     list(m = mClean, nipals_pca(mClean, k), indexRemaining = indexRem)
#   }, error = function(e) {
#     if(grepl(pattern = paste0("replacement has length zero"), x = e)) {
#       cleanParam <- cleanParam + (1 - cleanParam) / 2
#       message(paste("Too many NA in the data. Cleaning with maxNAs at", 
#                     cleanParam))
#       # pass the original m so indexRemaining is valid for the user's matrix
#       autoNipals(m, k, cleanParam)
#     } else { stop(e) }
#   })
# }



#### threshold ####
#' Apply threshold to a score or loading matrix
#'
#' Returns a binary matrix of the same size as \code{m} where all elements over
#' the threshold are 1.
#'
#' If th is a vector, the first element of th will be used as threshold for the
#' first col/row in m, etc.
setGeneric("threshold", signature = c("m", "th"), function(m, th, ...) {standardGeneric("threshold")})

#' @export
setMethod("threshold", c(m = "matrix", th = "numeric"), function(m, MARGIN = 2, th) {
  if(length(th) == 1) {
    mat <- matrix(TRUE, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
    mat[m < th] <- FALSE
    mat
  }
  else {
    if(MARGIN == 1) {
      if(length(th) != nrow(m)) { stop("Length of th must equal nrow(m).") }
      mat <- do.call(rbind, lapply(seq_len(nrow(m)), function(row) {
        m[row, ] > th[row]
      }))
    } else {
      if(length(th) != ncol(m)) { stop("Length of th must equal ncol(m)") }
      mat <- do.call(cbind, lapply(seq_len(ncol(m)), function(col) {
        m[, col] > th[col]
      }))
    }
    colnames(mat) <- colnames(m)
    rownames(mat) <- rownames(m)
    mat
  }
}
)

setMethod("threshold", c(m = "matrix", th = "matrix"), function(m, MARGIN = 2, th) {
  threshold(m, MARGIN, as.numeric(th))
}
)

validateBiclustNames <- function(biclustNames, bcs) {
  if (length(biclustNames) > 0) {
    validNames <- names(bcs)
    if (any(!biclustNames %in% validNames)) {
      warning(
        paste0(
          "Requested bicluster labels ",
          paste(setdiff(biclustNames, validNames), sep = ", "),
          " do not match any names in the \"pred\" slot of ",
          "the named BiclusterStrategy. Call \"?names\"",
          "to see how to view bicluster names."
        )
      )
      biclustNames <- intersect(biclustNames, validNames)
    }
    biclustNames
  } else {
    biclustNames
  }
}


validateK <- function(k) {
  if (!is.atomic(k)) {
    warning(
      "Argument \"k\" contains multiple elements. Attempting to use the
      first \"k\" provided."
    )
    k <- k[1]
  }
  
  if (!is.numeric(k)) {
    warning("Attempting to coerce \"k\" to numeric. Please inspect output
            carefully!")
    k <- as.numeric(k)
    if (is.na(k)) {
      stop("Couldn't coerce \"k\" to numeric.")
    }
  }
  
  if (!all(is.wholenumber(k))) {
    warning("Rounding \"k\" to whole number. Please inspect output carefully!")
    k <- round(k)
  }
  
  if (!(k > 0)) {
    stop("Argument \"k\" must be positive.")
  }
  
  k
}

validateKM <- function(k, m = NULL, method) {
  k <- validateK(k)
  if(method == "als-nmf" || "method" == "svd-pca" || method == "nipals-pca") {
    if (k >= min(nrow(m), ncol(m))) {
      warning(paste("Initializing k to the size of the smaller matrix",
                    "dimension."))
      return(min(nrow(m), ncol(m)))
    }
  }
  if(method == "plaid" || method == "spectral") k
  k
}

validateKs <- function(ks) {
  unlist(lapply(ks, function(k)
    validateK(k)))
}

validatePhenoNames <- function(phenoNames, bce) {
  if (length(phenoNames) > 0) {
    # If provided, phenoLabels must be in phenoData labels.
    validNames <- colnames(Biobase::phenoData(bce))
    if (any(!phenoNames %in% validNames)) {
      warning(
        paste0(
          "Requested phenoLabels ",
          paste(setdiff(phenoNames, validNames), sep = ", "),
          " are not in the phenoData slot of object x."
        )
      )
      phenoNames <- intersect(phenoNames, validNames)
    }
    phenoNames
  } else {
    phenoNames
  }
}


validateStratName <- function(stratName, bce) {
  validNames <- names(bce)
  if (length(stratName) > 1) {
    stop("More than one BiclusterStrategy name provided.")
  }
  if (length(stratName) == 0) {
    stop("BiclusterStrategy name was expected but not provided.")
  }
  if (!stratName %in% validNames) {
    stop(
      paste0(
        "Argument \"strategy\" is not the name of any",
        "BiclusterStrategy object in the provided ",
        "BiclusterExperiment."
      )
    )
  }
}

