#' Import sequencing data, sample metadata, and gene list
#'
#' @param expr_path File path to gene expression matrix (genes as rownames, samples as column names, cell fill as gene expression level). This gene expression matrix should be normalised (e.g., variance-stabilised).
#' @param sample_meta_path File path to sample metadata (one row per sample. column names of this table should include "sample", "treatment". within the sample column, the names should match the column names of the expression matrix)
#' @param gene_meta_path File path to gene metadata (one row per gene. column names of this table should include "gene", "category". within the gene column, the names should feature within the row names of the expression matrix)
#' @return A named list with elements: `expr`, `sample_meta`, and `gene_meta`.
#' @export
import_data <- function(expr_path, sample_meta_path, gene_meta_path) {
  expr <- read.csv(expr_path, row.names = 1, check.names = FALSE)
  sample_meta <- read.csv(sample_meta_path)
  gene_meta <- read.csv(gene_meta_path)

  if (!all(sample_meta$sample%in% colnames(expr))) {
    stop("Sample names in metadata do not match expression matrix column names.")
  }
  return(list(expr = expr, sample_meta = sample_meta, gene_meta = gene_meta))
}
