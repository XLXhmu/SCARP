# Declare variables used via non-standard evaluation (ggplot2/data.table).
utils::globalVariables(c("gene", "x", "y"))

#' Plot the cumulative enrichment curve of a single cell
#'
#' Visualizes the cumulative enrichment curve for a specified analysis unit.
#'
#' @param mat A numeric expression matrix with genes as rows and analysis
#'   units, such as samples, cells, or spatial spots, as columns.
#' @param cell_id A character string specifying the identifier of a single cell.
#' @param risk_score A numeric risk score corresponding to the cell specified
#'   by \code{cell_id}, as calculated by the \code{\link{Scarp}} function.
#' @param topN An integer specifying the number of top-ranked genes used to
#'   calculate the cumulative enrichment score. Default: 3000.
#' @param target_genes A character vector specifying the genes to display.
#'   Default: \code{NULL}.
#'
#' @returns A \code{ggplot} object displaying the cumulative enrichment curve.
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
#'
#' p <- cumulativeCurve(
#'   mat = SeuratObject::LayerData(object, layer = "counts"),
#'   cell_id = "AAACCTGTCACCCGAG-1_CRC14_tissue",
#'   risk_score = res[["AAACCTGTCACCCGAG-1_CRC14_tissue"]],
#'   topN = 3000,
#'   target_genes = c(
#'     "CD44", "CD24", "EPCAM", "PROM1", "SIRT4", "MET", "CXCR4",
#'     "SOX9", "SOX2", "NANOG", "ZEB1", "SNAI2", "HMGA1"
#'   )
#' )
#' print(p)
#' }
#'
#' @importFrom data.table data.table
#' @importFrom ggplot2 ggplot aes geom_line geom_area geom_hline geom_vline
#'   geom_point annotate theme_classic xlab ylab theme element_text
#'   element_line
#' @importFrom ggrepel geom_text_repel
#' @export
cumulativeCurve <- function(mat,
                            cell_id,
                            risk_score,
                            topN = 3000,
                            target_genes = NULL) {
  # ---- Input validation -----------------------------------------------------
  mat <- extract_expression_matrix(mat)
  gene_Diff <- load_internal_gene_diff()

  if (is.character(cell_id) && !(cell_id %in% colnames(mat))) {
    stop("`cell_id` not found in the columns of the expression data.")
  }
  if (!is.numeric(risk_score) || length(risk_score) != 1L) {
    stop("`risk_score` must be a single numeric value.")
  }
  topN <- as.integer(topN)
  if (length(topN) != 1L || topN < 1L) {
    stop("`topN` must be a single positive integer.")
  }

  genes <- intersect(rownames(mat), names(gene_Diff))
  if (length(genes) == 0) {
    stop("No genes in the expression data match the names of `gene_Diff`.")
  }

  mat <- mat[genes, , drop = FALSE]
  expr_vec <- as.numeric(as.matrix(mat[, cell_id, drop = FALSE]))
  names(expr_vec) <- genes

  order_idx <- order(expr_vec, decreasing = TRUE)
  topN <- min(topN, length(order_idx))
  topN_genes <- genes[order_idx[seq_len(topN)]]

  cumulative_scores <- cumsum(gene_Diff[topN_genes])

  plot_data_tmp <- data.table::data.table(
    x = seq_len(length(cumulative_scores)),
    y = cumulative_scores,
    gene = names(cumulative_scores)
  )

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

  # Colour scheme: red for positive risk, blue for negative, grey for zero.
  col_line <- "#808080"
  col_area <- "#c0c0c0"
  if (risk_score > 0) {
    col_line <- "#c12c1f"
    col_area <- "#ee7959"
  } else if (risk_score < 0) {
    col_line <- "#2e59a7"
    col_area <- "#88abda"
  }

  target_df <- NULL
  if (!is.null(target_genes)) {
    target_df <- plot_data_tmp[gene %in% target_genes]
  }

  plot1 <- ggplot2::ggplot(plot_data_tmp, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(linewidth = 0.235, color = col_line) +
    ggplot2::geom_area(fill = col_area, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "#999999", linewidth = 0.235) +
    ggplot2::geom_vline(xintercept = max_pos, linetype = "dashed",
                        color = "#999999", linewidth = 0.235) +
    ggplot2::geom_point(
      data = plot_data_tmp[x == max_pos],
      ggplot2::aes(x = x, y = y),
      color = "#8A0000", size = 1.5, shape = 16, stroke = 1
    )

  if (!is.null(target_df) && nrow(target_df) > 0) {
    plot1 <- plot1 +
      ggplot2::geom_point(
        data = target_df, ggplot2::aes(x = x, y = y),
        color = "#c12c1f", size = 1, shape = 16, stroke = 0.5
      ) +
      ggrepel::geom_text_repel(
        data = target_df,
        ggplot2::aes(x = x, y = y, label = gene),
        size = 2,
        color = "black",
        segment.color = "black",
        segment.size = 0.3,
        box.padding = 0.3,
        point.padding = 0.2,
        min.segment.length = 0.1,
        force = 1,
        force_pull = 1,
        max.overlaps = Inf,
        direction = "both",
        nudge_x = 0.1,
        nudge_y = 0.1,
        show.legend = FALSE
      )
  }

  plot1 <- plot1 +
    ggplot2::annotate(
      "text",
      x = max(plot_data_tmp$x) * 0.9,
      y = max(plot_data_tmp$y) * 0.9,
      label = paste0(
        "Cellid: ", cell_id, "\n",
        "Risk score: ", round(risk_score, 1), "\n",
        "Max deviation: ", round(max_dev, 1), "\n",
        "Rank: ", max_pos, "\n"
      ),
      hjust = 1, vjust = 1, size = 2, color = "black"
    ) +
    ggplot2::theme_classic() +
    ggplot2::ylab("Cumulative score") +
    ggplot2::xlab("Gene rank") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 6, colour = "black"),
      axis.title.y = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.text.y = ggplot2::element_text(size = 6, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.235),
      axis.ticks = ggplot2::element_line(linewidth = 0.235, colour = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 8),
      legend.position = "none"
    )

  plot1
}
