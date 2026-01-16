#' Dimensional reduction
#'
#' Dimensional reduction of counts using 0-1 scaling, PCA, or LDA.
#'
#' For each biological category (e.g., group of genes related to a given biological process), reduces the dimensionality across features to a single value for each sample.
#'
#' If using method="scale", the function returns the average value of the 0-1 scaled features within each category for each sample
#'
#' If using method="pca", the function uses a Principal Component Analysis (PCA) on the counts matrix to derive the coordinates of the samples for each group. The top PC dimensions are retained based on a variance threshold (default: threshold = 0.5). The axis which maximises the variance in sample-averaged PCA coordinates for the user-defined focus (default: focus = "group") is defined. Each sample is projected onto this axis. The most extreme samples are set to have values of 0 and 1, with lower projection values corresponding to the side of the axis showing lower counts values.
#'
#' If using method="lda", the function uses a Principal Component Analysis (PCA) on the counts matrix to derive the coordinates of the samples for each group. The top PC dimensions are retained based on a variance threshold (default: threshold = 0.5). When more than 1 PC is needed to reach threshold, a Linear Discriminant Analysis (LDA) is performed to increase the separation between samples based on a user-defined lda_focus (default: lda_focus = "group"). The top LD dimensions are retained based on a variance lda_threshold (default: lda_threshold = 0.8). The axis which maximises the variance in sample-averaged LDA coordinates for the user-defined focus (default: focus = "group") is defined. Each sample is projected onto this axis. The most extreme samples are set to have values of 0 and 1, with lower projection values corresponding to the side of the axis showing lower counts values.
#'
#' @param data_input Counts data, sample information, and biological process list uploaded using import_data()
#' @param method Choice of scale, pca, or lda
#' @param pca_threshold when using method = "pca" or "lda", threshold of cumulative variance used to identify the number of PCs to keep,
#' @param lda_threshold when using method = "lda", threshold of cumulative variance used to identify the number of LDs to keep,
#' @param pca_scale parameter used inside of prcomp() function. Default = TRUE.
#' @param focus when using method = "pca" or "lda", grouping of samples used to determine the main axis of variance (default = "group")
#' @param lda_focus when using method = "lda", grouping of samples used to perform LDA analysis (default = "group")
#' @param correlation_method statistical test used to determine linear correlation between values extracted for each sample and mean scaled counts in the category (default = "spearman").
#' @param remove_effect_from use limma to remove the effect of one or a few factors on the counts matrix before analysis (default: none. Use with "treatment" or c("treatment1","treatment2)).

#' @return When using method = "scale", the function returns a list with:
#'   - projection: data frame with sample, normalised_distance (normalized avg expr), category, furthest_groups.
#'   - information: data frame with number of features per category.
#'   When using method = "pca", the function returns a list with:
#'   - projection: data frame with sample, biological category, furthest sample groups based on the top PC dimensions retained, distance of each sample along a multidimensional line linking the most extreme samples, 0-1 scaled distance with direction based on counts of most extreme groups.
#'   - information: data frame with category, number of features (genes/proteins/others), method used ("PCA"), number of PCs kept, sum of variance kept, correlation between values extracted for samples in that biological category and their counts.
#'   - dimred_information: data frame with information needed for plot_dimensions(): biological category, number of PCs retained, percentage of variance described by PC1 and PC2, centroid, and slope of main axis of variance.
#'   When using method = "lda", the function returns a list with:
#'   - projection: data frame with sample, biological category, furthest sample groups based on the top PC dimensions retained, distance of each sample along a multidimensional line linking the most extreme samples, 0-1 scaled distance with direction based on counts of most extreme groups.
#'   - information: data frame with category, number of features (genes/proteins/others), method used ("PCA" or "LDA"), number of PCs and LDs kept, sum of variance kept, correlation between values extracted for samples in that biological category and their counts.
#'   - dimred_information: data frame with information needed for plot_dimensions(): biological category, number of PCs retained, percentage of variance described by PC1, PC2, LD1, and LD2, centroid, and slope of main axis of variance.
#' @export
dim_reduction <- function(data_input, method = c("scale","pca","lda"),pca_threshold=0.5,lda_threshold=0.8, pca_scale = TRUE,focus="group",lda_focus="group",correlation_method="spearman",remove_effect_from=NULL) {
  if (length(method) != 1) {
    stop("Please select exactly one method: 'scale', 'pca', 'lda'. Note that pca and lda will only run for categories with at least 5 features (e.g., genes).")
  }
  if (!requireNamespace("dplyr",quietly=TRUE)) stop("Package 'dplyr' is required.",call.=FALSE)

  counts = data_input$counts
  feature_meta = data_input$feature_meta
  sample_meta = data_input$sample_meta

  # Ensure rownames exist
  if (is.null(rownames(sample_meta))) {
    stop("rownames of the sample metadata must correspond to the column names from the counts matrix. \n")
  }

  # Reorder metadata to match columns of counts_mat
  if (!all(colnames(counts) %in% rownames(sample_meta))) {
    warning("Some samples from the counts matrix (column names) were not found in the sample metadata. \n")
  }

  sample_meta = sample_meta[which(sample_meta$sample%in%colnames(counts)), ]

  # Warn if the order was changed
  if (!all(rownames(sample_meta) == colnames(counts))) {
    warning("Sample metadata was reordered to match the column order from the counts matrix. \n")
  }

  sample_meta = sample_meta[match(colnames(counts),sample_meta$sample),]
  data_input2 = data_input

  if(length(remove_effect_from)>0){
      design = model.matrix(as.formula(paste("~ ",paste(remove_effect_from,collapse="+"))),data = sample_meta)
      fit = limma::lmFit(counts, design)
      counts = limma::residuals.MArrayLM(fit, counts)
      print(paste0("The counts data was adjusted based on the following: ",paste("~ ",paste(remove_effect_from,collapse="+"))))
      data_input2$counts = counts
  }
  method <- match.arg(method, choices = c("pca", "scale","lda"))


  if (method == "pca") {
    if(length(which(colnames(sample_meta)%in%focus))==0){
      stop("focus (default: group) column not found in sample information table. Make sure to have a column named group or use argument focus.")
    }
    res <- pca_method(data_input2,pca_threshold,focus,correlation_method,pca_scale)
  } else if (method == "scale") {
    res <- scale_method(data_input2)
  } else if (method =="lda"){
    if(length(which(colnames(sample_meta)%in%focus))==0){
      stop("focus (default: group) column not found in sample information table. Make sure to have a column named group or use argument focus.")
    }
    res <- lda_method(data_input2, pca_threshold, lda_threshold,focus,lda_focus,correlation_method,pca_scale)

  }

  return(res)
}

