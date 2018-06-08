genBiclusters <- function(sizeX, sizeY, biclusterRows, biclusterCols, noise = 0) {
  if (length(biclusterRows) != length(biclusterCols)) {
    stop("biclusterRows and biclusterCols must be the same length")
  }
  
  res <- matrix(0, nrow = sizeX, ncol = sizeY)
  
  invisible(mapply(function(rowRange, colRange) { res[rowRange, colRange] <<- 1 },
         rowRange = biclusterRows, colRange = biclusterCols))
  
  # multiply every column in res by a random integer in the set [1,5]
  k <- runif(n = sizeY, min = 1, max = 5)
  res <- res * rep(k, rep.int(nrow(res), length(k)))
  
  res <- res + noise * matrix(rnorm(n = sizeX * sizeY), nrow = sizeX)
  res <- res - min(res)
  
  #plot
  old.par <- par(no.readonly=T) 
  par(mar=c(0, 0, 0, 0))
  image(res, useRaster=TRUE, axes=FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), c(min(res), round(max(res), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  res
}

genSimData <- function() {
  res <- genBiclusters(50, 250, rowIds, colIds, 0.25)
  res
}