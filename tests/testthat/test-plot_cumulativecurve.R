test_that("cumulativeCurve returns a ggplot for matrix input", {
  mat <- load_sce_counts(100)
  res <- Scarp(mat, topN = 50)
  p <- cumulativeCurve(mat, cell_id = colnames(mat)[1],
                       risk_score = res[1], topN = 50)

  expect_s3_class(p, "ggplot")
})

test_that("cumulativeCurve works with a Seurat object", {
  seurat <- make_seurat(60)
  res <- Scarp(seurat, topN = 50)
  p <- cumulativeCurve(seurat, cell_id = colnames(seurat)[1],
                       risk_score = res[1], topN = 50)

  expect_s3_class(p, "ggplot")
})

test_that("cumulativeCurve handles a zero risk score", {
  mat <- load_sce_counts(100)
  p <- cumulativeCurve(mat, cell_id = colnames(mat)[1],
                       risk_score = 0, topN = 50)

  expect_s3_class(p, "ggplot")
})

test_that("cumulativeCurve handles NULL target_genes", {
  mat <- load_sce_counts(100)
  p <- cumulativeCurve(mat, cell_id = colnames(mat)[1],
                       risk_score = -1.5, topN = 50, target_genes = NULL)

  expect_s3_class(p, "ggplot")
})

test_that("cumulativeCurve errors on an unknown cell_id", {
  mat <- load_sce_counts(100)

  expect_error(
    cumulativeCurve(mat, cell_id = "nope", risk_score = 1,
                    topN = 50),
    "`cell_id` not found"
  )
})

test_that("cumulativeCurve validates risk_score", {
  mat <- load_sce_counts(100)

  expect_error(
    cumulativeCurve(mat, cell_id = colnames(mat)[1], risk_score = c(1, 2),
                    topN = 50),
    "`risk_score` must be a single numeric value"
  )
})

test_that("cumulativeCurve validates topN", {
  mat <- load_sce_counts(100)

  expect_error(
    cumulativeCurve(mat, cell_id = colnames(mat)[1], risk_score = 1,
                    topN = 0),
    "`topN` must be a single positive integer"
  )
})
