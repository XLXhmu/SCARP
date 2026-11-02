test_that("reference retains supplied NR-versus-R values", {
  data("gene_Diff", package = "SCARP")
  expect_length(gene_Diff, 19157L)
  expect_true(all(is.finite(gene_Diff)))
  expect_false(anyDuplicated(names(gene_Diff)) > 0L)
  expect_true(all(diff(gene_Diff) <= 0))
  expect_equal(unname(gene_Diff["VAX1"]), 4.762449, tolerance = 1e-6)
})

test_that("normalization applies the specified position and direction weights", {
  ans <- SCARP:::calculate_enrichment_score(c(2, 5, -1),
                                           random_mean = 1, random_sd = 2,
                                           nGenes = 10)
  expect_equal(ans$normalized_score, ((5-1)/2) * (1-2/10)^2 * (5/6))
  expect_equal(ans$max_position, 2L)
})

test_that("background excludes the lower half and preserves sample size", {
  sig <- setNames(c(1, 2, 100, 1000), letters[1:4])
  ans <- SCARP:::calculate_deviation(sig, nperm = 20,
                                    gene_list = letters[1:4], topN = 2)
  expect_equal(ans$random_mean, 3)
  expect_error(SCARP:::calculate_deviation(sig, nperm = 20,
                gene_list = letters[1:4], topN = 3), "top-50%")
})

test_that("scoring is reproducible and dense and sparse paths agree", {
  mat <- load_sce_counts(6)
  set.seed(123); sparse <- Scarp(mat, topN = 50)
  set.seed(123); dense <- Scarp(as.matrix(mat), topN = 50)
  expect_equal(sparse, dense)
  expect_true(all(is.finite(sparse)))
  mat[, 1] <- 0
  expect_equal(unname(Scarp(mat, topN = 50)[1]), 0)
})

test_that("invalid topN inputs receive an informative error", {
  mat <- load_sce_counts(2)
  for (value in list(NA_real_, Inf, 1.5, numeric(), c(1,2), "50")) {
    expect_error(Scarp(mat, topN = value), "single positive integer")
  }
  expect_error(Scarp(mat, topN = nrow(mat)), "top-50%")
})
