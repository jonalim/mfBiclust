context("BiclusterStrategy")
set.seed(12345)

# takes a long time, how can we pare this down?
test_that("the BiclusterStrategy constructor is functional", {
  m <- matrix(c(25, 49, 20, 72, 38, 44, 43, 39, 26, 61, 1, 7, 53, 64, 34),
              nrow = 3, ncol = 5)
  bcs_snmf <- BiclusterStrategy(m, k = 2, bicluster = "pca",
                                scoreThresh = "otsu",
                                loadingThresh = "otsu")
  expect_true(validObject(bcs_snmf))
})

test_that("snmf handles weird matrices", {
  set.seed(1)
  m <- matrix(1:16, nrow = 4)
  expect_warning({bcs <- BiclusterStrategy(m, k = 3, bicluster = "snmf/l",
                                           scoreThresh = "otsu",
                                           loadingThresh = "otsu")},
                 "Sparse NMF failed, switching to PCA.")
})

test_that("thresholding is accurate", {
  load("../testdata/ref_threshold.rda")
  res <- generateThresholdMatrix("otsu", l[[1]], l[[2]])
  expect_equal(res, l[[3]])
})