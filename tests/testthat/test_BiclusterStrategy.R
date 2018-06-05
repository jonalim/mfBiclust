context("BiclusterStrategy")
library(NMF)
set.seed(12345)
r <- 3
c <- 5
k <- 2
m <- matrix(c(25, 49, 20, 72, 38, 44, 43, 39, 26, 61, 1, 7, 53, 64, 34),
            nrow = r, ncol = c)
bcs_snmf <- BiclusterStrategy(m, k = k, bicluster = "snmf/l",
                         scoreThresh = "otsu",
                         loadingThresh = "otsu")

bcs_pca <- BiclusterStrategy(m, k = k, bicluster = "pca",
                         scoreThresh = "otsu",
                         loadingThresh = "otsu")

test_that("BiclusterStrategy catches illegal input", {
	expect_error(BiclusterStrategy(m, 3))
	})

test_that("matrix factors have expected dimensions", {
  expect_equal(dim(score(bcs_snmf)), c(r, k))
  expect_equal(dim(loading(bcs_snmf)), c(k, c))
  expect_equal(dim(score(bcs_pca)), c(r, k))
  expect_equal(dim(loading(bcs_pca)), c(k, c))
})

test_that("matrix factorization is accurate", {
  load("../testdata/ref_pca.rda")
  load("../testdata/ref_snmf.rda")
  expect_equal(score(bcs_snmf), ref_snmf@factors@fit@W)
  expect_equal(loading(bcs_snmf), ref_snmf@factors@fit@H)
  expect_equal(bcs_snmf@scoreThresh, ref_snmf@scoreThresh)
  expect_equal(score(bcs_pca), ref_pca@factors@fit@W)
  expect_equal(loading(bcs_pca), ref_pca@factors@fit@H)
  expect_equal(bcs_pca@scoreThresh, ref_pca@scoreThresh)
})
          
test_that("snmf handles weird matrices", {
  m <- matrix(1:16, nrow = 4)
  bcs <- BiclusterStrategy(m, k = 3, bicluster = "snmf/l",
                           scoreThresh = "otsu",
                           loadingThresh = "otsu")
})