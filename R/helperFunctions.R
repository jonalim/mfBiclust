#' @include generics.R
NULL

#' Convert logical matrices to two lists of rows and columns
#'
#' Given two logical matrices \eqn{A_{m,k}} and \eqn{B_{k,n}} showing which rows
#' and columns are in bicluster k, returns two \code{\link{list}s}. The first
#' lists the rows in each bicluster. The second lists the columns in each
#' bicluster.
#'
#' @param rowxBicluster an m x k logical matrix
#' @param biclusterxCol a k x n logical matrix
#'
#' @return A \code{\link{list}} containing two \code{list}s
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' biclusterLists <- biclusterMatrix2List(clusteredFeatures(bcs),
#'  clusteredSamples(bcs))
#' biclusterLists[[1]] # Rows in each bicluster
#' biclusterLists[[2]] # Cols in each bicluster
#' @export
biclusterMatrix2List <- function(rowxBicluster, biclusterxCol) {
  biclusterRows <- lapply(seq_len(ncol(rowxBicluster)),
                          function(k) which(rowxBicluster[, k]))
  biclusterCols <- lapply(seq_len(nrow(biclusterxCol)),
                          function(k) which(biclusterxCol[k, ]))
  list(biclusterRows, biclusterCols)
}

#' Convert row and column numbers to a matrix representation
#'
#' This function converts the output from
#' \code{\link[biclust:bicluster]{biclusternumber}} into the membership matrix
#' format. The first returned membership matrix shows which rows are in which
#' bicluster. The second returned membership matrix shows which columns are in
#' which bicluster.
#'
#' @param biclusterNumber A list describing the rows and columns in each
#'   bicluster
#' @param m The numeric matrix in which biclusters were computed
#' @param k Restricts the results to the first \code{k} biclusters
#'
#' @return a \code{\link{list}} of two matrices
#'
#' @export
#'
#' @examples
#' require(biclust)
#' data <- matrix(rnorm(400), 20, 20)
#' data[12:16, 8:12] <- rnorm(25, 3, 0.3)
#' set.seed(1)
#' biclusters <- biclust(data, BCPlaid(), back.fit = 2, shuffle = 3,
#'  fit.model = ~m + a + b, iter.startup = 5, iter.layer = 30, verbose = TRUE)
#' biclusterNumbers <- biclusternumber(biclusters)
#' biclusterNumbers2membershipMatrices(biclusterNumbers, data, 1)
biclusterNumbers2membershipMatrices <- function(biclusterNumber, m, k) {
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
  if (s == "als-nmf") { return("ALS-NMF") }
  if (s == "nipals-pca") { return("NIPALS-PCA") }
  if (s == "svd-pca") { return("SVD-PCA") }
  else { switch(s, snmf = return("SNMF"), pca = return("PCA"),
                otsu = return("Otsu"),
                return(paste0(toupper(substring(s, 1,1)), substring(s, 2)))) }
})

#### Clean ####
setMethod("clean", c(object = "matrix"), function(object, cleanParam,
                                                  dimsRemain) {
  if(!(cleanParam <= 1 && cleanParam >= 0)) {
    stop("Arg \"cleanParam\" must be in the range of 0 to 1.")
  }
  maxNaPerRow <- round((1 - cleanParam) * ncol(object))
  maxNaPerCol <- round((1 - cleanParam) * nrow(object))

  # Both 0s and NAs can foul up NIPALS
  goodRows <- apply(object, MARGIN = 1, function(row)
    (sum(is.na(row)) + sum(row == 0, na.rm = TRUE)) < maxNaPerRow)
  goodCols <- apply(object, MARGIN = 2, function(col)
    (sum(is.na(col)) + sum(col == 0, na.rm = TRUE)) < maxNaPerCol)

  object <- object[goodRows, goodCols]

  if(dimsRemain) {
    list(obj = object, dimsRemain = list(goodRows, goodCols))
  } else { object }

})

# Use a string to set the random seed
#
# Returns the output of .Random.seed. Call
# \code{on.exit(assign(".Random.seed", oldSeed, envir=globalenv()), add = TRUE)
# } to restore the random seed.
duplicable <- function(str) {
  if (!exists(".Random.seed", mode="numeric")) sample(NA)
  oldSeed <- .Random.seed
  newSeed <- strtoi(str, 35L)
  if(!is.na(newSeed)) set.seed(newSeed) else {
    set.seed(12345)
    warning("Argument could not be encoded as numeric. Seed set to 12345")
  }
  oldSeed
}

#' Filter biclusters by overlap and quantity
#'
#' First all biclusters encompassing the whole matrix are removed. Next,
#' biclusters are selected in descending order of size, skipping any biclusters
#' that overlap excessively with already- selected biclusters. The remaining
#' \eqn{k^\prime} biclusters are returned in membership-matrix format, along
#' with a matrix giving the union of the retained biclusters and a logical
#' vector telling which biclusters were retained. See
#' \code{\link{clusteredSamples}()} for a description of the membership matrix
#' format.
#'
#' @param rowxBicluster a row-bicluster membership matrix
#' @param biclusterxCol a bicluster-column membership matrix
#' @param max the maximum number of biclusters to keep
#' @param overlap the maximum portion of a bicluster that may overlap with
#'   another bicluster. Passed directly to \code{\link{overlap}()}
#'
#' @return A list of the format:
#'  \describe{
#'    \item{RowxBicluster}{A row-bicluster membership matrix
#'    \eqn{A_{m,k^\prime}}, where \eqn{k^\prime} is the number of retained
#'    biclusters}
#'    \item{BiclusterxCol}{A bicluster-column membership matrix
#'    \eqn{B_{k^\prime,n}}}
#'    \item{biclustered}{A binary numeric matrix giving the union of all
#'    retained biclusters}
#'    \item{chosen}{A logical vector of length \eqn{k}}
#'  }
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 4, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' filtered <- filter.biclust(clusteredFeatures(bcs), clusteredSamples(bcs),
#'   overlap = 0.25)
#' filtered$chosen
#' # The second bicluster was excluded because .31 of its matrix elements
#' # overlapped with the first bicluster.
#'
#' @export
filter.biclust <- function(rowxBicluster, biclusterxCol, max = NULL,
                           overlap = 0.25) {

    # Could create a version of this function that takes a BiclusterStrategy as
    # input and returns a modified BiclusterStrategy. Would fit with a redesign
    # of BiclusterStrategy that contains a record of all functions called to
    # modify it.
  if(ncol(rowxBicluster) != nrow(biclusterxCol)) { # validate
    stop(paste0("RowxBicluster must have the same number of columns as",
                "BiclusterxCol"))
  }

  k <- ncol(rowxBicluster)
  if(k == 0) {
    chosen = rep(FALSE, ncol(rowxBicluster))
  } else if(k == 1 ) {
    # no filtering needed
    chosen = rep(TRUE, ncol(rowxBicluster))
  } else {
    # Create lists of rows and columns contained in biclusters
    rc <- biclusterMatrix2List(rowxBicluster, biclusterxCol)
    biclusterRows <- rc[[1]]
    biclusterCols <- rc[[2]]

    chosen <- rep(FALSE, each = k)
    pool <- rep(TRUE, each = k)
    sizes <- sizes(biclusterRows, biclusterCols)
    names(sizes) <- seq_len(k)

    # exclude empty biclusters and whole-dataset biclusters
    pool[sizes == 0 | sizes == nrow(rowxBicluster) * ncol(biclusterxCol)] <- FALSE
    # For each bicluster, calculate its overlap with all other biclusters
    overlaps <- overlap(biclusterRows, biclusterCols, TRUE)

        # Start with no biclusters
    # Evaluates to TRUE even if argument max was missing
    while(all(sum(chosen) < max) && sum(pool) > 0) {
        # Choose biclusters in order of size, decreasing
      chooseMe <- as.numeric(names(which.max(sizes[which(pool)])))
      chosen[chooseMe] <- TRUE

      # Remove from the pool any biclusters heavily overlapping with the chosen
      # bicluster. This prevents any two biclusters overlapping more than overlap
      pool[overlaps[, chooseMe] > overlap] <- FALSE
      pool[chooseMe] <- FALSE
    }
  }
  # We've chosen enough when we've chosen as many as possible without violating
  # the overlap limit, OR when we've chosen max biclusters.

  # Determine the union of all biclustered matrix elements
  biclustered <- union.biclust(rowxBicluster, biclusterxCol, chosen)

  list(RowxBicluster = rowxBicluster[, chosen, drop = FALSE],
       BiclusterxCol = biclusterxCol[chosen, , drop = FALSE],
       biclustered = biclustered, chosen = chosen)
}

is.wholenumber <-
  function(x, tol = sqrt(.Machine$double.eps)) {
    abs(x - round(x)) < tol
  }

logicalMatrix2Numeric <- function(m) {
  matrix(as.numeric(m), dim(m), dimnames = dimnames(m))
}

#' Calculate overlaps between every pair of biclusters.
#'
#' The overlap between biclusters \eqn{i} and \eqn{j} is defined as the number
#' of matrix elements shared by the two biclusters divided by the size of
#' bicluster \eqn{i}.
#'
#' @param BiclusterRows a list, each element containing row indices in the
#'  bicluster
#' @param BiclusterCols a list, each element containing column indices in the
#'  bicluster
#' @param matrix if true, returns a matrix \eqn{A} where element \eqn{A_{i,j}}
#'   is the overlap between biclusters \eqn{i} and \eqn{j}. Otherwise returns a
#'   list where the \eqn{i}th element is a vector giving the overlap between
#'   bicluster \eqn{i} and each bicluster, including itself. In both instances,
#'   the order of biclusters is preserved from the input.
#' @export
#'
#' @return If \code{matrix}, a numeric matrix. Else, a list of numeric vectors.
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' biclusterLists <- biclusterMatrix2List(clusteredFeatures(bcs),
#'  clusteredSamples(bcs))
#' biclusterRows <- biclusterLists[[1]]
#' biclusterCols <- biclusterLists[[2]]
#' overlap(biclusterRows, biclusterCols)
#' overlap(biclusterRows, biclusterCols, matrix = TRUE)
overlap <- function(BiclusterRows, BiclusterCols, matrix = FALSE) {
  sizes <- vapply(seq_along(BiclusterRows), FUN.VALUE = numeric(1), FUN = function(biclus) {
    length(BiclusterRows[[biclus]]) * length(BiclusterCols[[biclus]])
  })

  overlaps <- lapply(seq_along(BiclusterRows), function(biclus1) {
    vapply(seq_along(BiclusterRows), FUN.VALUE = numeric(1), FUN = function(biclus2) {
      intersection1 <- base::intersect(BiclusterRows[[biclus1]], BiclusterRows[[biclus2]])
      intersection1 <- if(length(intersection1) > 0) length(intersection1)
      else 0
      intersection2 <- base::intersect(BiclusterCols[[biclus1]], BiclusterCols[[biclus2]])
      intersection2 <- if(length(intersection2) > 0) length(intersection2)
      else 0
      intersection <- intersection1 * intersection2
      overlap <- intersection / sizes[biclus1]
    })
  })
  if(matrix) {
    return(as.matrix(do.call(rbind, overlaps)))
  } else {
    return(overlaps)
  }
}

pseudovalues <- function(m) {
  if(any(m < 0)) {
    message(paste("Converting to pseudovalues (x + abs(min(x))) just for",
                  "this BiclusterStrategy because",
                  "negative values are not allowed."))
    m <- m + abs(min(m))
  }
  return(m)
}



#' Calculate bicluster sizes
#'
#' Calculates the number of matrix elements in each bicluster. Biclusters must
#' be described in the same format as the output from
#' \code{\link{biclusterMatrix2List}}.
#'
#' @param biclusterRows a list of numeric vectors giving the row indices of each
#'  bicluster
#' @param biclusterCols a list of numeric vectors giving the column indices of
#'  each bicluster
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 2, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#' biclustersAsLists <- biclusterMatrix2List(
#'  rowxBicluster = clusteredFeatures(bcs),
#'  biclusterxCol = clusteredSamples(bcs))
#' sizes(biclusterRows = biclustersAsLists[[1]],
#'       biclusterCols = biclustersAsLists[[2]])
sizes <- function(biclusterRows, biclusterCols) {
  size <- vapply(seq_along(biclusterRows), FUN.VALUE = numeric(1), FUN = function(biclus) {
    length(biclusterRows[[biclus]]) * length(biclusterCols[[biclus]])
  })
  return(size)
}

#' Create a binary matrix describing a union of biclusters
#'
#' Given feature-membership matrix \eqn{A_{m,k}} and sample-membership matrix
#' \eqn{B_{k,n}} describing a set of biclusters, returns a binary matrix where
#' any element that is a member of any bicluster is 1.
#'
#' @param RowxBiclust a boolean row-membership matrix
#' @param BiclustxCol a boolean column-membership matrix
#' @param choose If provided, each element of \code{choose} corresponds to a
#'   bicluster (a column of \code{RowxBiclust} and a row of \code{BiclustxCol}.
#'   FALSE elements will exclude biclusters from the calculation.
#'
#' @return numeric matrix
#' @export
#'
#' @examples
#' bce <- BiclusterExperiment(yeast_benchmark[[1]])
#' bce <- addStrat(bce, k = 4, method = "als-nmf")
#' bcs <- getStrat(bce, 1)
#'
#' # Combine all biclusters
#' union.biclust(clusteredFeatures(bcs), clusteredSamples(bcs))
#'
#' # Combine only the first and second biclusters
#' union.biclust(clusteredFeatures(bcs), clusteredSamples(bcs),
#'  c(TRUE, TRUE, FALSE, FALSE))
union.biclust <- function(RowxBiclust, BiclustxCol,
                          choose = rep(TRUE, ncol(RowxBiclust))) {
  biclustered <- matrix(0, nrow = nrow(RowxBiclust),
                        ncol = ncol(BiclustxCol))
  lapply(which(choose), function(bicluster) {
    biclustered[RowxBiclust[ , bicluster], BiclustxCol[bicluster, ]] <<- 1
  })
  biclustered
}

validateBiclustNames <- function(biclustNames, bcs) {
  if (length(biclustNames) > 0) {
    validNames <- bcNames(bcs, allBc = TRUE)
    if (any(!biclustNames %in% validNames)) {
      warning(
        paste0(
          "Requested bicluster labels ",
          paste(setdiff(biclustNames, validNames), sep = ", "),
          " do not match any of bcNames(bcs)")
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

# Ensure k is valid for the selected matrix and method
#
# Modify this when adding a new biclustering method
validateKM <- function(k, m = NULL, method) {
  k <- validateK(k)
  if(method == "als-nmf" || "method" == "svd-pca" || method == "nipals-pca" ||
     method == "snmf") {
    if (k > min(nrow(m), ncol(m))) {
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

