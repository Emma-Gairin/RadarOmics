#' Import sequencing data, sample metadata, and molecular feature list
#'
#' @param counts_path File path to counts matrix (molecular features as rownames, samples as column names, cell fill as counts level). This counts matrix should be normalised (e.g., variance-stabilised).
#' @param sample_meta_path File path to sample metadata (one row per sample. column names of this table should include "sample", "treatment". within the sample column, the names should match the column names of the counts matrix)
#' @param feature_meta_path File path to molecular feature (e.g., gene) metadata (one row per feature. column names of this table should include "feature", "category". within the molecular feature column, the names should feature within the row names of the counts matrix)
#' @return A named list with elements: `counts`, `sample_meta`, and `feature_meta`.
#' @export
import_data <- function(counts_path, sample_meta_path, feature_meta_path) {
  counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  sample_meta <- read.csv(sample_meta_path)
  feature_meta <- read.csv(feature_meta_path)
  print(paste0(length(intersect(rownames(counts),feature_meta$feature))," features were retrieved in the row names of the counts matrix."))
  if(length(intersect(rownames(counts),feature_meta$feature))==0){
    stop("No molecular feature from the biological information match the rownames of the counts matrix. \n Make sure that that the first column of your counts matrix CSV file is named 'feature' and corresponds to the molecular feature names/IDs.")
  }
  print(paste0(length(intersect(colnames(counts),sample_meta$sample))," samples were retrieved in the column names of the counts matrix."))

  if (!all(sample_meta$sample%in% colnames(counts))) {
    warning("Sample names in metadata do not match counts matrix column names. We only kept the ones in common")
    sample_meta = sample_meta[which(sample_meta$sample%in%colnames(counts)),]
  }
  return(list(counts = counts, sample_meta = sample_meta, feature_meta = feature_meta))
}
