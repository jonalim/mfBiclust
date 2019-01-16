set.seed(12345)
# No background, sixteen non-overlapping 5x5 biclusters. Elements in
# biclusters have value 5; all other elements are 0.
fiveSquareConstant <- genSimData(n = 16, biclusterConstant = 5,
 clusterHeight = 5, clusterWidth = 5, shuffle = FALSE)

# No background, two 5x5 biclusters causing row-shift effects
fiveSquareRowShift <- genSimData(n = 2, rowBase = 1, rowShift = 1,
 clusterHeight = 5, clusterWidth = 5, shuffle = FALSE)

# Three 5x5 biclusters, where the first and second biclusters overlap and the
# second and third biclusters overlap. Both overlap regions are 4 rows by 2
# columns.
fiveSquareOverlapped <- genSimData(n = 3, rowBase = 1, rowShift = 1, overlapRows = 4,
 overlapCols = 2, clusterHeight = 5, clusterWidth = 5, shuffle = FALSE)

# One 10x10 plaid bicluster
plaid <- genSimData(n = 3, biclusterShift = 1, rowShift = 1, colShift = 1,
 clusterHeight = 10, clusterWidth = 10, shuffle = FALSE)

simdata <- list(fiveSquareConstant, fiveSquareRowShift,
                fiveSquareOverlapped, plaid)
names(simdata) <- c("fiveSquareConstant", "fiveSquareRowShift",
                    "fiveSquareOverlapped", "plaid")
usethis::use_data(simdata, overwrite = TRUE)
