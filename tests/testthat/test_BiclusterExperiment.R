context("BiclusterExperiment")

input <- matrix(c(25, 49, 20, 72, 38, 44, 43, 39, 26, 61, 1, 7, 53, 64, 34),
                nrow = 3, ncol = 5)
bce <- BiclusterExperiment(m = input)

test_that("BiclusterExperiment constructor works", {
  expect_true(validObject(bce))
})