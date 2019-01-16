#' Simulate an omics dataset with biclusters
#'
#' Most models of omics datasets can be built using various combinations of
#' effects in this function. Effects are applied in this order:
#' \enumerate{
#' \item Constant background value
#' \item Base row values (overrides constant background value)
#' \item Gaussian noise
#' \item Constant value across all biclusters (overrides all previous effects)
#' \item Additive bicluster-specific values
#' \item Additive bicluster-and-column-specific values
#' \item Multiplicative bicluster-and-column-specific values
#' \item Additive bicluster-and-row-specific values
#' \item Multiplicative bicluster-and-row-specific values
#' \item Shuffling of rows and columns
#' }
#'
#' Besides constant values across all biclusters (\code{biclusterConstant}),
#' all effects are sampled from a normal distribution centered on 0 and with
#' standard deviation given by the respective parameter.
#'
#' @return A numeric matrix containing \code{n} biclusters created by the given
#'   parameters.
#'
#' @param n the number of biclusters
#' @param clusterHeight the number of rows in each bicluster
#' @param clusterWidth the number of columns in each bicluster
#' @param dimx the number of rows in the dataset
#' @param dimy the number of columns in the dataset
#' @param overlapRows the number of rows that biclusters should share
#' @param overlapCols the number of columns that biclusters should share
#' @param bgConst initial value of the entire matrix
#' @param rowBase standard deviaton of base row values
#' @param bgNorm standard deviation of Gaussian noise
#' @param biclusterConstant if not NULL, sets the value of every element in
#'   every bicluster, after creating the background
#' @param biclusterShift standard deviation of bicluster-specific additive
#'  values
#' @param colShift standard deviation of bicluster-and-column-specific shifts
#' @param colScale standard deviation of bicluster-and-column-specific scaling
#'   coefficients
#' @param rowShift standard deviation of bicluster-and-row-specific shifts
#' @param rowScale standard deviation of bicluster-and-row-specific scaling
#'   coefficients
#' @param shuffle if FALSE, biclusters will be deterministically arranged
#'   along the diagonal of the matrix
#' @param file an optional string giving a prefix for .png and .csv files
#'   where the generated matrix should be written
#'
#' @export
#'
#' @examples
#' # No background, sixteen non-overlapping 5x5 biclusters. Elements in
#' # biclusters have value 5; all other elements are 0.
#' mat <- genSimData(n = 16, biclusterConstant = 5,
#'  clusterHeight = 5, clusterWidth = 5, shuffle = FALSE)
#'
#' # No background, two 5x5 biclusters causing row-shift effects
#' mat <- genSimData(n = 2, rowBase = 1, rowShift = 1,
#'  clusterHeight = 5, clusterWidth = 5, shuffle = FALSE)
#' matrixHeatmap(mat)
#'
#' # Three 5x5 biclusters, where the first and second biclusters overlap and the
#' # second and third biclusters overlap. Both overlap regions are 4 rows by 2
#' # columns.
#' mat <- genSimData(n = 3, rowBase = 1, rowShift = 1, overlapRows = 4,
#'  overlapCols = 2, clusterHeight = 5, clusterWidth = 5, shuffle = FALSE)
#' matrixHeatmap(mat)
#'
#' # One 10x10 plaid bicluster
#' mat <- genSimData(n = 3, biclusterShift = 1, rowShift = 1, colShift = 1,
#'  clusterHeight = 10, clusterWidth = 10, shuffle = FALSE)
#'  matrixHeatmap(mat)
genSimData <- function(n = 1, clusterHeight = 20, clusterWidth = 20,
                       dimx = 80, dimy = 80,
                       overlapRows = 0, overlapCols = 0,
                       bgConst = 0, bgNorm = 0, biclusterConstant = NULL,
                       biclusterShift = 0, rowBase = 0, rowShift = 0,
                       rowScale = NULL, colShift = 0, colScale = NULL,
                       shuffle = TRUE, file = "") {
  if(!(
    (length(clusterHeight) == 1 || length(clusterHeight) == n) &&
    (length(clusterWidth) == 1 || length(clusterWidth) == n) &&
    all(overlapRows <= clusterHeight) &&
    all(overlapCols <= clusterWidth) &&
    clusterHeight + (clusterHeight - overlapRows) * (n - 1) <= dimx &&
    clusterWidth + (clusterWidth - overlapCols) * (n - 1) <= dimy
    )) {
    stop(paste("Incompatible dimensions. Arguments must satisfy:\n",
               "(length(clusterHeight) == 1 || length(clusterHeight) == n) &&
                 (length(clusterWidth) == 1 || length(clusterWidth) == n) &&
                 all(overlapRows <= clusterHeight) &&
                 all(overlapCols <= clusterWidth) &&
                 clusterHeight + (clusterHeight - overlapRows) * (n - 1)",
               "<= dimx &&
                 clusterWidth + (clusterWidth - overlapCols) * (n - 1)",
               "<= dimy\n"))
  }

  xStarts <- (clusterHeight - overlapRows) * (seq_len(n) - 1) + 1
  xStops <- xStarts + clusterWidth - 1
  yStarts <- (clusterWidth - overlapCols) * (seq_len(n) - 1) + 1
  yStops <- yStarts + clusterHeight - 1

  xCoords <- mapply(`:`, xStarts, xStops, SIMPLIFY = FALSE)
  yCoords <- mapply(`:`, yStarts, yStops, SIMPLIFY = FALSE)

  res <- genSimDataHelper(sizeX = dimx, sizeY = dimy,
                          biclusterConstant = biclusterConstant,
                          biclusterShift = biclusterShift,
                          biclusterRows = xCoords,
                          biclusterCols = yCoords,
                          rowBase = rowBase, rowShift = rowShift, rowScale = rowScale,
                          colShift = colShift,
                          colScale = colScale,
                          bgConst = bgConst, bgNorm = bgNorm)

  if(shuffle) {
    res <- res[sample(seq_len(dimx), size = dimx, replace = FALSE),
               sample(seq_len(dimy), size = dimy, replace = FALSE)]
  }
  if (nchar(file) > 0) {
    grDevices::png(paste0(file, ".png"))
    write.csv(res, row.names = FALSE, file = paste0(file, ".csv"))
    }

  #plot
  old.par <- par(no.readonly = TRUE)
  par(mar = c(0, 0, 0, 0))
  image(t(apply(res, 2, rev)), useRaster = TRUE, axes = FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(x = "topright", legend = c(min(res), round(max(res), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  if (nchar(file) > 0) { grDevices::dev.off() }
  return(res)
}

genSimDataHelper <- function(sizeX, sizeY, biclusterRows, biclusterCols,
                             biclusterConstant,
                             biclusterShift,
                             rowBase, rowShift, rowScale, colScale, colShift,
                             bgConst, bgNorm) {
  if (length(biclusterRows) != length(biclusterCols)) {
    stop("biclusterRows and biclusterCols must be the same length")
  }

  res <- matrix(bgConst, nrow = sizeX, ncol = sizeY)

  # Background is row-stripes
  invisible(lapply(seq_len(nrow(res)), function(row) {
    res[row, ] <<- rnorm(1, 0, rowBase)
  }))

  # add Gaussian noise
  res <- res + bgNorm * matrix(rnorm(n = sizeX * sizeY), nrow = sizeX)

  # Every bicluster set to the same value
  if(!is.null(biclusterConstant)) {
      invisible(
          mapply(function(rowRange, colRange) {
              res[rowRange, colRange] <<- biclusterConstant
          },
          rowRange = biclusterRows, colRange = biclusterCols))
  }

  # make biclusters additive if desired
  invisible(mapply(function(rowRange, colRange) {
    res[rowRange, colRange] <<- res[rowRange, colRange] + rnorm(1, 0, biclusterShift)

    # add to every bicluster column a random number in the set [1, colShift],
    # normally distributed
    k <- rnorm(n = length(colRange), mean = 0, sd = colShift)
    res[rowRange, colRange] <<- res[rowRange, colRange] + rep(k, rep.int(length(rowRange), length(k)))

    # multiply every column in the bicluster by a random number in the set
    # [1, colScale], normally distributed
    if(!is.null(colScale)) {
      k <- rnorm(n = length(colRange), mean = 0, sd = colScale)
      res[rowRange, colRange] <<- res[rowRange, colRange] *
        rep(k, rep.int(length(rowRange), length(k)))
    }

    # same for rows
    k <- rnorm(n = length(rowRange), mean = 0, sd = rowShift)
    res[rowRange, colRange] <<- res[rowRange, colRange] + rep(k, times = length(colRange))

    if(!is.null(rowScale)) {
      k <- rnorm(n = length(rowRange), mean = 0, sd = rowScale)
      res[rowRange, colRange] <<- res[rowRange, colRange] * rep(k, times = length(colRange))
    }
  },
  rowRange = biclusterRows, colRange = biclusterCols))

  res
}
