context("BiclusterExperiment")

oldSeed <- mfBiclust:::duplicable("test")
input <- matrix(runif(n = 36), nrow = 6)
input[1:2, 1:2] <- 20
bce <- BiclusterExperiment(m = input)

singular <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 4)
bce.singular <- BiclusterExperiment(singular)

# How I generated the reference instances of BiclusterExperiment
# ref_als_nmf <- addStrat(bce, k = 1, method = "als-nmf")
# ref_svd_pca <- addStrat(bce, k = 1, method = "svd-pca")
# ref_snmf <- addStrat(bce, k = 1, method = "snmf")
# addStrat(bce.singular, k = 2, method = "snmf", silent = TRUE)
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
# don't compare BiclusterExperiment@.__classVersion__
test_that("ALF-NMF works", {
    expect_equivalent(addStrat(bce, k = 1, method = "als-nmf"), bce.als_nmf)
})
test_that("SVD-PCA works", {
    expect_equivalent(addStrat(bce, k = 1, method = "svd-pca"), bce.svd_pca)
})

test.snmf <- addStrat(bce, k = 1, method = "snmf")
test_that("SNMF is accurate (sample biclustering)", {
    # snmf records its runtime in the NMFfit object in BiclusterStrategy@fit
    # we only care about the hard-bicluster matrices; all other portions
    # of the pipeline are tested in other tests
    expect_equivalent(clusteredSamples(getStrat(test.snmf, 1)),
                      clusteredSamples(getStrat(bce.snmf, 1)))
})
test_that("SNMF is accurate (feature biclustering)", {
    expect_equivalent(clusteredFeatures(getStrat(test.snmf, 1)),
                      clusteredFeatures(getStrat(bce.snmf, 1)))
})
# test_that("SNMF handles singular matrices correctly", {
#     expect_warning(addStrat(bce.singular, k = 2, method = "snmf",
#                             silent = TRUE),
#                    regexp = "snmf failed, switching to PCA")
# })

test_that("NIPALS-PCA works", {
    expect_equivalent(addStrat(bce, k = 1, method = "nipals-pca",
                               silent = TRUE),
                      bce.nipals_pca)
})
test_that("Plaid works", {
    expect_equivalent(addStrat(bce, k = 1, method = "plaid", silent = TRUE),
                      bce.plaid)
})
test_that("Spectral works", {
    expect_warning(test.spectral <- addStrat(bce, k = 1, method = "spectral",
                                             silent = TRUE))
    expect_equivalent(test.spectral, bce.spectral)
})

assign(".Random.seed", oldSeed, envir=globalenv())
