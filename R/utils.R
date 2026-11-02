# Internal utility functions for SCARP -----------------------------------------
#
# These functions are not exported and are used internally by Scarp()
# and cumulativeCurve().


#' Coerce the input to a gene-by-cell expression matrix
#'
#' Accepts either a dense \code{matrix}, a sparse \code{Matrix}, or a Seurat
#' object (from which the RNA count matrix is extracted).
#'
#' @param mat Expression matrix, sparse matrix, or Seurat object.
#' @returns A matrix-like object with row names (genes).
#' @noRd
extract_expression_matrix <- function(mat) {
  if (inherits(mat, "Seurat")) {
    mat <- SeuratObject::LayerData(mat, layer = "counts")
  }
  if (is.null(rownames(mat)) || is.null(ncol(mat))) {
    stop("`mat` must be a gene-by-cell expression matrix or a Seurat object.")
  }
  if (nrow(mat) < 1L || ncol(mat) < 1L) {
    stop("`mat` must have at least one gene and one cell.")
  }
  if (anyDuplicated(rownames(mat))) {
    stop("Gene row names must be unique.")
  }
  if (is.null(colnames(mat)) || anyDuplicated(colnames(mat))) {
    stop("Cell column names must be present and unique.")
  }
  values <- if (inherits(mat, "sparseMatrix")) mat@x else as.vector(as.matrix(mat))
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("Expression values must be numeric and finite.")
  }
  mat
}


#' Load the bundled risk signature
#'
#' @returns A named numeric vector of embedded gene-level risk statistics.
#' @noRd
load_internal_gene_diff <- function() {
  local_env <- environment()
  utils::data("gene_Diff", package = "SCARP", envir = local_env)
  gene_Diff <- get("gene_Diff", envir = local_env, inherits = FALSE)
  if (!is.numeric(gene_Diff) || is.null(names(gene_Diff))) {
    stop("The bundled `gene_Diff` dataset must be a named numeric vector.")
  }
  gene_Diff
}


#' Build the permutation-based null distribution
#'
#' Randomly draws \code{topN} genes (from a trimmed "null universe") for
#' \code{nperm} permutations, computes the maximum signed deviation of the
#' cumulative sum of \code{gene_Diff}, and returns the mean and standard
#' deviation of the null distribution.
#'
#' @param gene_Diff Named numeric vector of per-gene differential values.
#' @param nperm Number of permutations.
#' @param gene_list Character vector of gene names (ranked by mean expression)
#'   used to build the null universe.
#' @param topN Number of genes drawn in each permutation.
#' @returns A list with elements \code{random_mean} and \code{random_sd}.
#' @importFrom stats sd
#' @noRd
calculate_deviation <- function(gene_Diff, nperm = 3000, gene_list, topN) {
  n_genes <- length(gene_list)

  # Draw exactly topN genes from the upper half of matched genes.
  trim_size <- floor(n_genes / 2)
  if (topN > trim_size || trim_size < 1L) {
    stop("`topN` exceeds the top-50% background pool.")
  }
  gene_list <- gene_list[seq_len(trim_size)]

  ES_null <- numeric(nperm)
  for (i in seq_len(nperm)) {
    sorted_indices <- sort(sample.int(trim_size, topN))
    topN_genes <- gene_list[sorted_indices]
    cumulative_scores <- cumsum(gene_Diff[topN_genes])
    max_positive <- max(cumulative_scores)
    max_negative <- min(cumulative_scores)
    max_dev <- if (abs(max_positive) >= abs(max_negative)) {
      max_positive
    } else {
      max_negative
    }
    ES_null[i] <- max_dev
  }

  random_sd <- stats::sd(ES_null)
  if (is.na(random_sd) || random_sd == 0) {
    random_sd <- .Machine$double.eps
  }

  list(random_mean = mean(ES_null), random_sd = random_sd)
}


#' Compute a normalised enrichment score from a cumulative curve
#'
#' Internal helper that turns a cumulative-sum vector into a normalised
#' enrichment score by comparing its maximum signed deviation against a
#' permutation-based null model.
#'
#' @param cumulative_scores Numeric vector of cumulative scores.
#' @param random_mean Mean of the null distribution (default 0).
#' @param random_sd Standard deviation of the null distribution (default 1).
#' @param nGenes Total number of genes used for position weighting.
#' @returns A list with the normalised score and summary statistics of the
#'   cumulative curve.
#' @noRd
calculate_enrichment_score <- function(cumulative_scores,
                                       random_mean = 0,
                                       random_sd = 1,
                                       nGenes) {
  max_positive <- max(cumulative_scores)
  max_negative <- min(cumulative_scores)

  if (abs(max_positive) >= abs(max_negative)) {
    direction <- "positive"
    max_dev <- max_positive
    max_pos <- which.max(cumulative_scores)
  } else {
    direction <- "negative"
    max_dev <- max_negative
    max_pos <- which.min(cumulative_scores)
  }

  position_weight <- (1 - (max_pos / nGenes))^2

  total_deviation <- abs(max_positive) + abs(max_negative)
  if (total_deviation == 0) {
    direction_strength <- 0
  } else {
    direction_strength <- if (direction == "positive") {
      abs(max_positive) / total_deviation
    } else {
      abs(max_negative) / total_deviation
    }
  }

  relative_dev <- (max_dev - random_mean) / random_sd
  normalized_score <- relative_dev * position_weight * direction_strength

  n <- length(cumulative_scores)
  positive_area <- sum(cumulative_scores > 0) / n
  negative_area <- sum(cumulative_scores < 0) / n

  list(
    normalized_score = normalized_score,
    direction = direction,
    max_deviation = max_dev,
    max_position = max_pos,
    position_weight = position_weight,
    direction_strength = direction_strength,
    positive_area = positive_area,
    negative_area = negative_area,
    cumulative_scores = cumulative_scores
  )
}
