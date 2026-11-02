#' Single-Cell Chemoresistance Assessment via Rank-Based Projection
#'
#' We developed a rank-based cumulative scoring algorithm, termed SCARP, to quantify chemotherapy resistance at single-cell resolution in colorectal cancer (CRC). Gene-specific log fold-change (logFC) values were derived from differential expression analysis between chemotherapy non-response (NR) and response (R) samples in a CRC chemotherapy-response dataset and were incorporated into SCARP as a reference chemoresistance signature. Under the NR-versus-R contrast, positive logFC values represent genes preferentially expressed in non-response samples, whereas negative logFC values represent genes preferentially expressed in response samples.
#' 
#' For each cell, genes were ranked according to their normalized expression levels, and the corresponding gene-specific logFC values were cumulatively summed along the ranked gene list. The cumulative deviation with the largest absolute magnitude was defined as the observed score (MaxScore). To account for expression-dependent background variation, TopN genes were randomly sampled without replacement from the top 50% of matched genes ranked by mean expression for 3,000 iterations. Sampled genes retained their mean-expression rank order. The resulting empirical distribution was used to estimate the background mean (Random_mean) and standard deviation (Random_sd).
#' 
#' The observed score was standardized against this background and further adjusted using a position weight (Pw), which emphasizes deviations arising from highly ranked genes, and a direction strength term (Ds), which captures the dominance of the non-response- or response-associated cumulative trajectory. The final single-cell chemoresistance score was calculated as (MaxScore - Random_mean) / Random_sd * Pw * Ds. Positive scores indicate a more positive signed deviation relative to the empirical background, associated with the non-response signature; negative scores indicate a more negative deviation relative to that background. Because the score is background-centered, its sign need not equal the sign of the unstandardized cumulative deviation. It is a transcriptional signature score, not a probability of clinical chemoresistance.
#'
#' @param mat A numeric expression matrix with genes as rows and analysis
#'   units, such as samples, cells, or spatial spots, as columns.
#' @param topN An integer specifying the number of top-ranked genes used to
#'   calculate the cumulative enrichment score. Default: 3000.
#'
#' @returns A named numeric vector of normalized enrichment scores, with names
#'   corresponding to the analysis unit identifiers.
#'
#' @examples
#' \dontrun{
#' # Load the example data
#' data(object)
#'
#' # Compute risk scores
#' res <- Scarp(
#'   mat = SeuratObject::LayerData(object, layer = "counts")
#' )
#' }
#'
#' @importFrom stats runif sd
#' @importFrom SeuratObject LayerData
#' @export
Scarp <- function(mat, topN = 3000) {
  mat <- extract_expression_matrix(mat)
  gene_Diff <- load_internal_gene_diff()

  if (!is.numeric(topN) || length(topN) != 1L || is.na(topN) ||
      !is.finite(topN) || topN < 1 || topN != floor(topN)) {
    stop("`topN` must be a single positive integer.")
  }

  topN <- as.integer(topN)

  genes <- intersect(rownames(mat), names(gene_Diff))
  nGenes <- length(genes)
  if (nGenes == 0) {
    stop("No genes in the expression data match the names of `gene_Diff`.")
  }

  if (topN > floor(nGenes / 2)) {
    stop("`topN` must not exceed the top-50% background pool (",
         floor(nGenes / 2), " matched genes). Use a smaller `topN`.")
  }

  mat <- mat[genes, , drop = FALSE]
  gene_Diff <- gene_Diff[genes]

  # ---- Build the permutation-based null model -------------------------------
  # Row means via sparse-safe matrix multiplication (rowMeans() does not
  # dispatch reliably on sparse matrices inside the package namespace).
  gene_means <- as.vector(mat %*% rep(1, ncol(mat))) / ncol(mat)
  names(gene_means) <- rownames(mat)
  gene_list <- names(sort(gene_means, decreasing = TRUE))

  perm_estimates <- calculate_deviation(
    gene_Diff,
    nperm = 3000,
    gene_list = gene_list,
    topN = topN
  )

  # ---- Per-cell scoring (chunked) ------------------------------------------
  num_cells <- ncol(mat)
  AUCs <- numeric(num_cells)
  names(AUCs) <- colnames(mat)

  chunk_size <- 2000
  num_chunks <- ceiling(num_cells / chunk_size)

  message("Processing ", num_cells, " cells in ", num_chunks, " chunk(s)...")

  for (c in seq_len(num_chunks)) {
    start_col <- (c - 1L) * chunk_size + 1L
    end_col <- min(c * chunk_size, num_cells)

    mat_chunk <- as.matrix(mat[, start_col:end_col, drop = FALSE])

    for (i in seq_len(ncol(mat_chunk))) {
      expr_vec <- mat_chunk[, i]

      non_zero_idx <- which(expr_vec != 0)

      if (length(non_zero_idx) == 0) {
        AUCs[start_col + i - 1L] <- 0
        next
      }

      expr_non_zero <- expr_vec[non_zero_idx]

      # Random tie-breaking; set.seed() makes repeated calls reproducible.
      noise <- stats::runif(length(expr_non_zero))
      order_idx <- order(-expr_non_zero, noise)

      if (length(order_idx) > topN) {
        order_idx <- order_idx[seq_len(topN)]
      }

      original_row_idx <- non_zero_idx[order_idx]

      cumulative_scores <- cumsum(gene_Diff[original_row_idx])

      result <- calculate_enrichment_score(
        cumulative_scores,
        random_mean = perm_estimates$random_mean,
        random_sd = perm_estimates$random_sd,
        nGenes = nGenes
      )

      AUCs[start_col + i - 1L] <- result$normalized_score
    }

    message(sprintf("  chunk %d/%d done.", c, num_chunks))
  }

  AUCs
}
