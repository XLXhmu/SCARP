test_that("Scarp returns a named numeric score vector", {
  mat <- load_sce_counts(100)
  res <- Scarp(mat, topN = 50)

  expect_type(res, "double")
  expect_named(res, colnames(mat))
  expect_length(res, ncol(mat))
})

test_that("Scarp accepts a Seurat object", {
  seurat <- make_seurat(60)
  res <- Scarp(seurat, topN = 50)

  expect_length(res, ncol(seurat))
  expect_named(res, colnames(seurat))
})

test_that("Scarp works on a dense matrix too", {
  mat <- load_sce_counts(60)
  mat_dense <- as.matrix(mat)
  res <- Scarp(mat_dense, topN = 50)

  expect_length(res, ncol(mat_dense))
})

test_that("Scarp errors when no genes match the signature", {
  mat <- matrix(0, nrow = 5, ncol = 3,
                dimnames = list(paste0("fake", 1:5), paste0("c", 1:3)))

  expect_error(Scarp(mat),
               "No genes in the expression data")
})

test_that("Scarp validates topN", {
  mat <- load_sce_counts(50)

  expect_error(Scarp(mat, topN = 0),
               "`topN` must be a single positive integer")
})
