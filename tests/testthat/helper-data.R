# Shared test helpers: load SCARP's example data and build a small
# Seurat object.

load_gene_diff <- function() {
  data("gene_Diff", package = "SCARP")
  gene_Diff
}

load_sce_counts <- local({
  cached <- NULL
  function(n_cells = 100) {
    if (is.null(cached)) {
      e <- new.env()
      data("object", package = "SCARP", envir = e)
      cached <<- SeuratObject::LayerData(e$object, layer = "counts")
    }
    cached[, seq_len(min(n_cells, ncol(cached))), drop = FALSE]
  }
})

make_seurat <- function(n_cells = 100) {
  counts <- load_sce_counts(n_cells)
  suppressWarnings(
    SeuratObject::CreateSeuratObject(counts = counts, project = "test")
  )
}
