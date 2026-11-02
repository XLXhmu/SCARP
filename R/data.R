#' CRC chemotherapy-response reference signature
#'
#' Named logFC values supplied from a Nature Medicine CRC dataset,
#' calculated by DESeq2 for non-response (NR) versus response (R).
#' Positive values are associated with non-response; negative values
#' are associated with response. Values are sorted in descending order.
#' @format A named numeric vector of length 19157; names are gene symbols.
#' @source User-provided file: Nature_medicine NRvsR DEG logFC.rda.
#'   Article DOI, accession, and DESeq2 settings have not been supplied.
#' @examples
#' data(gene_Diff)
#' head(gene_Diff)
"gene_Diff"

#' Pretreatment CRC single-cell example data
#'
#' User-provided Zhang Zemin single-cell data from two patients,
#' sampled from tumor tissue before treatment.
#' @format A Seurat object with 33538 genes and 10419 cells.
#' @source User-provided obj_zhang_sub.rds; article accession not supplied.
#' @examples
#' data(object)
#' dim(object)
"object"
