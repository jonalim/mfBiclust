context("BiclusterExperiment")

# input <- matrix(c(25, 49, 20, 72, 38, 44, 43, 39, 26, 61, 1, 7, 53, 64, 34),
                # nrow = 3, ncol = 5)
oldSeed <- duplicable("test")
input <- matrix(runif(n = 36), nrow = 6)
input[1:2, 1:2] <- 20
bce <- BiclusterExperiment(m = input)

# How I generated the reference instances of BiclusterExperiment
# ref_als_nmf <- addStrat(bce, k = 1, method = "als-nmf")
# ref_svd_pca <- addStrat(bce, k = 1, method = "svd-pca")
# ref_snmf <- addStrat(bce, k = 1, method = "snmf")
# ref_nipals_pca <- addStrat(bce, k = 1, method = "nipals-pca")
# ref_plaid <- addStrat(bce, k = 1, method = "plaid")
# ref_spectral <- addStrat(bce, k = 1, method = "spectral")
# ref_bces <- list(bce.als_nmf = ref_als_nmf,
#                  bce.svd_pca = ref_svd_pca,
#                  bce.snmf = ref_snmf,
#                  bce.nipals_pca = ref_nipals_pca,
#                  bce.plaid = ref_plaid, bce.spectral = ref_spectral)
# save(ref_bces, file = "../testdata/ref_bces.rda")

load(file = "../testdata/ref_bces.rda")

list2env(x = ref_bces, envir = environment())

test_that("BiclusterExperiment constructor works", {
  expect_true(validObject(bce))
})
test_that("ALF-NMF works", {
  expect_equal(addStrat(bce, k = 1, method = "als-nmf"), bce.als_nmf)
})
test_that("SVD-PCA works", {
  expect_equal(addStrat(bce, k = 1, method = "svd-pca"), bce.svd_pca)
})
test_that("SNMF works", {
  # snmf records its runtime in the NMFfit object in BiclusterStrategy@fit
  # we only care about the hard-bicluster matrices; all other portions
  # of the pipeline are tested in other tests
  test.snmf <- addStrat(bce, k = 1, method = "snmf")
  expect_equal(clusteredSamples(getStrat(test.snmf, 1)), 
               clusteredSamples(getStrat(bce.snmf, 1)))
  expect_equal(clusteredFeatures(getStrat(test.snmf, 1)),
               clusteredFeatures(getStrat(bce.snmf, 1)))
})
test_that("NIPALS-PCA works", {
  expect_equal(addStrat(bce, k = 1, method = "nipals-pca"), bce.nipals_pca)
})
test_that("Plaid works", {
  expect_equal(addStrat(bce, k = 1, method = "plaid"), bce.plaid)
})
test_that("Spectral works", {
  expect_warning(test.spectral <- addStrat(bce, k = 1, method = "spectral"))
  expect_equal(test.spectral, bce.spectral)
})

assign(".Random.seed", oldSeed, envir=globalenv())
