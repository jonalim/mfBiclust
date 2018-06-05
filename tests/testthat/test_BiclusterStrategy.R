context("BiclusterStrategy")

r <- 3
c <- 5
k <- 2
m <- matrix(c(25, 49, 20, 72, 38, 44, 43, 39, 26, 61, 1, 7, 53, 64, 34),
            nrow = r, ncol = c)
bcs <- BiclusterStrategy(m, k = k, bicluster = "snmf/l",
                         scoreThresh = "otsu",
                         loadingThresh = "otsu")

test_that("BiclusterStrategy catches illegal input", {
	expect_error(BiclusterStrategy(m, 3))
	})

test_that("the constructor maintains data integrity", {
  expect_equal(bcs@biclustAlgo, "snmf/l")
  expect_equal(bcs@scoreThreshAlgo, "otsu")
  expect_equal(bcs@loadingThreshAlgo, "otsu")
})

test_that("matrix factors have expected dimensions", {
  expect_equal(dim(score(bcs)), c(r, k))
  expect_equal(dim(loading(bcs)), c(k, c))
})
          
test_that("snmf handles weird matrices", {
  m <- matrix(1:16, nrow = 4)
  bcs <- BiclusterStrategy(m, k = 3, bicluster = "snmf/l",
                           scoreThresh = "otsu",
                           loadingThresh = "otsu")
})