# Run from the SCARP package directory. Original user-provided files are retained.
input_dir <- if (dir.exists("data/pre")) "data/pre" else "data"
signature_file <- file.path(input_dir, "Nature_medicine NRvsR DEG logFC.rda")
object_file <- file.path(input_dir, "obj_zhang_sub.rds")
if (!all(file.exists(c(signature_file, object_file)))) {
  stop("Original input files are required only for data regeneration. ",
       "Place them in data/pre/ first. Installation uses the bundled .rda datasets.")
}
e <- new.env(parent = emptyenv())
load(signature_file, envir = e)
gene_Diff <- e$gene_Diff
stopifnot(is.numeric(gene_Diff), !is.null(names(gene_Diff)),
          all(is.finite(gene_Diff)), !anyDuplicated(names(gene_Diff)),
          all(nzchar(names(gene_Diff))), all(diff(gene_Diff) <= 0))
save(gene_Diff, file = "data/gene_Diff.rda", compress = "xz", version = 3)
suppressPackageStartupMessages(library(SeuratObject))
object <- readRDS(object_file)
stopifnot(inherits(object, "Seurat"))
save(object, file = "data/object.rda", compress = "xz", version = 3)
cat("Signature genes:", length(gene_Diff), "| example dimensions:", dim(object), "\n")
