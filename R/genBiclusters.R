#' Simulate a matrix with n biclusters
#'
#' Optional noise parameter
#'
#' @export
genSimData <- function(n, overlapped = FALSE, file = "",
                       striped = c("cols", ""), noise = 0.01, dynamSize = FALSE,
                       dimx = 50, dimy = dimx) {
  
  striped <- match.arg(striped)
  
  if(dynamSize) {
    # scale up in size if more biclusters
    x <- 5 * n
    y <- 5 * n
  } else {
    # cram all biclusters into the same space. scale dimensions manually
    x <- dimx
    y <- dimy
  }
  xCoords <- sample(1:x, 2 * n, replace = FALSE)
  yCoords <- sample(1:y, 2 * n, replace = FALSE)
  if(!overlapped) {
    if(runif(1) > 0.5) { xCoords <- sort(xCoords) }
    else {yCoords <- sort(yCoords)}
  }
  xCoords <- lapply(seq_len(n), function(i) {
    xCoords[i*2-1]:xCoords[i*2]
  })
  yCoords <- lapply(seq_len(n), function(i) {
    yCoords[i*2-1]:yCoords[i*2]
  })
  
  res <- genSimDataHelper(sizeX = x, sizeY = y, biclusterRows = xCoords,
                       biclusterCols = yCoords, striped = striped,
                       noise = noise)
  
  if (nchar(file) > 0) { 
    png(paste0(file, ".png"))
    write.csv(res, row.names = FALSE, file = paste0(file, ".csv"))}
  
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
                             noise = 0, striped = c("", "cols")) {
  if (length(biclusterRows) != length(biclusterCols)) {
    stop("biclusterRows and biclusterCols must be the same length")
  }
  striped <- match.arg(striped)
  
  res <- matrix(0, nrow = sizeX, ncol = sizeY)
  
  invisible(mapply(function(rowRange, colRange) {
    biclusterVal <- rep(1, each = length(colRange))
    if(striped == "cols") {
      # stripes will have a uniform distribution; range of 2
      colstripes <- runif(length(biclusterVal), 0, 2)
      biclusterVal <- biclusterVal + colstripes
    }
    mapply(function(col, value) {
      # overlapping biclusters will be additive
      res[rowRange, col] <<- res[rowRange, col] + value
      },
      col = colRange, value = biclusterVal)
  },
  rowRange = biclusterRows, colRange = biclusterCols))
  
  # multiply every column in res by a random integer in the set [1,5]
  k <- ceiling(runif(n = sizeY, min = 0, max = 5))
  res <- res * rep(k, rep.int(nrow(res), length(k)))
  
  res <- res + noise * matrix(rnorm(n = sizeX * sizeY), nrow = sizeX)
  res <- res - min(res)
  
  old.par <- par(no.readonly=T) 
  par(mar=c(0, 0, 0, 0))
  image(t(apply(res, 2, rev)), useRaster=TRUE, axes=FALSE, col = RColorBrewer::brewer.pal(9, "BuPu"))
  legend(grconvertX(0.5, "device"), grconvertY(1, "device"), c(min(res), round(max(res), digits = 3)),
         fill = RColorBrewer::brewer.pal(9, "BuPu")[c(1, 9)])
  par(old.par)
  
  res
}
