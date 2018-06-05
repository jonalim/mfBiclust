context("Annotation parsing")
bce <- mfbc(matrix(sample(1:10, size = 16, replace = TRUE), nrow = 3), nogui = TRUE)
# createAnnots <- function(x, names, strategy = "", phenoLabels = c(), biclustLabels = c()) {

test_that("Annotation parsing creates an annotation dataframe", {
  expect_equal(createAnnots(bce, names = 
  
})