#' Simulate a matrix with n biclusters
#'
#' Optional noise parameter
#'
#' @export
genSimData <- function(n = 1, clusterHeight = 20, clusterWidth = 20,
                       dimx = 80, dimy = 80, 
                       overlapRows = 0, overlapCols = 0,
                       biclusterConstant = NULL, biclusterShift = 0, 
                       rowShift = 0, rowScale = NULL,
                       colShift = 0, colScale = NULL,
                       bgConst = 0, bgNorm = 0, bgUnif = 0, 
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
                          rowShift = rowShift, rowScale = rowScale,
                          colShift = colShift,
                          colScale = colScale,
                          bgConst = bgConst, bgNorm = bgNorm, bgUnif = bgUnif)
  
  if(shuffle) {
    res <- res[sample(1:dimx, size = dimx, replace = FALSE),
               sample(1:dimy, size = dimy, replace = FALSE)]
  }
  if (nchar(file) > 0) { 
    png(paste0(file, ".png"))
    write.csv(res, row.names = FALSE, file = paste0(file, ".csv"))
    }
  
  #plot
  old.par <- par(no.readonly = T)
  par(mar = c(0, 0, 0, 0))
  image(t(apply(res, 2, rev)), useRaster = TRUE, axes = FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"),
         c(min(res), round(max(res), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  if (nchar(file) > 0) { dev.off() }
  return(res)
}

genSimDataHelper <- function(sizeX, sizeY, biclusterRows, biclusterCols, 
                             biclusterConstant,
                             biclusterShift,
                             rowShift, rowScale, colScale, colShift,
                             bgConst, bgNorm, bgUnif) {
  if (length(biclusterRows) != length(biclusterCols)) {
    stop("biclusterRows and biclusterCols must be the same length")
  }
  
  res <- matrix(bgConst, nrow = sizeX, ncol = sizeY)
  
  # add Gaussian noise
  res <- res + bgNorm * matrix(rnorm(n = sizeX * sizeY), nrow = sizeX)
  
  # add uniform background
  res <- res + matrix(runif(n = sizeX * sizeY, 0, bgUnif), nrow = sizeX)
  
  # Every bicluster set to the same value
  if(!is.null(biclusterConstant)) {
  invisible(mapply(function(rowRange, colRange) {
    res[rowRange, colRange] <<- biclusterConstant
  },
  rowRange = biclusterRows, colRange = biclusterCols))
  }
  
  # make biclusters additive if desired
  invisible(mapply(function(rowRange, colRange) {
    res[rowRange, colRange] <<- res[rowRange, colRange] + rnorm(1, 0, biclusterShift)

    # multiply every column in the bicluster by a random number in the set 
    # [1, colScale], normally distributed
    if(!is.null(colScale)) {
      k <- rnorm(n = length(colRange), mean = 0, sd = colScale)
      res[rowRange, colRange] <<- res[rowRange, colRange] * 
        rep(k, rep.int(length(rowRange), length(k)))
    }
    
    # add to every bicluster column a random number in the set [1, colShift],
    # normally distributed
    k <- rnorm(n = length(colRange), mean = 0, sd = colShift)
    res[rowRange, colRange] <<- res[rowRange, colRange] + rep(k, rep.int(length(rowRange), length(k)))

        # same for rows
    if(!is.null(rowScale)) {
      k <- rnorm(n = length(rowRange), mean = 0, sd = rowScale)
      res[rowRange, colRange] <<- res[rowRange, colRange] * rep(k, times = length(colRange))
    }
    
    k <- rnorm(n = length(rowRange), mean = 0, sd = rowShift)
    res[rowRange, colRange] <<- res[rowRange, colRange] + rep(k, times = length(colRange))
  },
  rowRange = biclusterRows, colRange = biclusterCols))
  
  old.par <- par(no.readonly=T) 
  par(mar=c(0, 0, 0, 0))
  image(t(apply(res, 2, rev)), useRaster=TRUE, axes=FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), c(min(res), round(max(res), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  
  res
}
