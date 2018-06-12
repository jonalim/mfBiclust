genBiclusters <- function(sizeX, sizeY, biclusterRows, biclusterCols, noise = 0) {
  if (length(biclusterRows) != length(biclusterCols)) {
    stop("biclusterRows and biclusterCols must be the same length")
  }
  
  res <- matrix(0, nrow = sizeX, ncol = sizeY)
  invisible(mapply(function(rowRange, colRange) { res[rowRange, colRange] <<- 1 },
         rowRange = biclusterRows, colRange = biclusterCols))
  
  # multiply every column in res by a random integer in the set [1,5]
  k <- ceiling(runif(n = sizeY, min = 0, max = 5))
  res <- res * rep(k, rep.int(nrow(res), length(k)))
  
  res <- res + noise * matrix(rnorm(n = sizeX * sizeY), nrow = sizeX)
  res <- res - min(res)
  
  old.par <- par(no.readonly=T) 
  par(mar=c(0, 0, 0, 0))
  image(res, useRaster=TRUE, axes=FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), c(min(res), round(max(res), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  
  res
}

#' @export
genSimData3 <- function(save = FALSE) {
  res <- genBiclusters(50, 250, list(10:25, 26:38, 41:45), list(25:159, 101:200, 201:250), 0.25)
  res
}

#' Simulate a matrix with n biclusters
#'
#' Optional noise parameter
#'
#' @export
genSimData <- function(n, noise = 0.01) {
  x <- round(rnorm(1, 20, 1) * n)
  y <- round(rnorm(1, 85, 5) * n)
  
  xCoords <- sort(sample(1:x, 2 * n, replace = FALSE)) 
  xCoords <- lapply(seq_len(n), function(i) {
    xCoords[i*2-1]:xCoords[i*2]
  })
  yCoords <- sort(sample(1:y, 2 * n, replace = FALSE)) 
  yCoords <- lapply(seq_len(n), function(i) {
    yCoords[i*2-1]:yCoords[i*2]
  })
  
  genBiclusters(sizeX = x, sizeY = y, biclusterRows = xCoords,biclusterCols = yCoords, noise = noise)
}
