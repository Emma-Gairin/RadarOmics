#' Dimensional reduction of expression counts using 0-1 scaled expression levels, PCA, or LDA
#'
#' For each category of variables (e.g., group of genes related to a given biological process), reduces the dimensionality across variables to a single value for each sample.
#'
#' If using method="scale", the function returns the average value of the 0-1 scaled variable within each category for each sample
#'
#' If using method="pca", the function uses a Principal Component Analysis (PCA) on the expression matrix to derive the coordinates of the samples for each group. The top PC dimensions are retained based on a variance threshold (default: threshold = 0.5). The axis which maximises the variance in sample-averaged PCA coordinates for the user-defined focus (default: focus = "group") is defined. Each sample is projected onto this axis. The most extreme samples are set to have values of 0 and 1, with lower projection values corresponding to the side of the axis showing lower expression values.
#'
#' If using method="lda", the function uses a Principal Component Analysis (PCA) on the expression matrix to derive the coordinates of the samples for each group. The top PC dimensions are retained based on a variance threshold (default: threshold = 0.5). When more than 1 PC is needed to reach threshold, a Linear Discriminant Analysis (LDA) is performed to increase the separation between samples based on a user-defined lda_focus (default: lda_focus = "group"). The top LD dimensions are retained based on a variance lda_threshold (default: lda_threshold = 0.8). The axis which maximises the variance in sample-averaged LDA coordinates for the user-defined focus (default: focus = "group") is defined. Each sample is projected onto this axis. The most extreme samples are set to have values of 0 and 1, with lower projection values corresponding to the side of the axis showing lower expression values.
#'
#' @param data_input Expression data, sample information, and biological process list uploaded using import_data()
#' @param method Choice of scale, pca, or lda
#' @param threshold when using method = "pca" or "lda", threshold of cumulative variance used to identify the number of PCs to keep,
#' @param lda_threshold when using method = "lda", threshold of cumulative variance used to identify the number of LDs to keep,
#' @param focus when using method = "pca" or "lda", grouping of samples used to determine the main axis of variance (default = "group")
#' @param lda_focus when using method = "lda", grouping of samples used to perform LDA analysis (default = "group")
#' @param correlation_method statistical test used to determine linear correlation between values extracted for each sample and mean scaled gene expression in the category (default = "spearman").
#' @param remove_effect_from use limma to remove the effect of one or a few variables on the expression matrix before analysis (default: none. Use with "treatment" or c("treatment1","treatment2)).

#' @return When using method = "scale", the function returns a list with:
#'   - projection: data frame with sample, normalised_distance (normalized avg expr), category, furthest_groups.
#'   - information: data frame with number of genes per category.
#'   When using method = "pca", the function returns a list with:
#'   - projection: data frame with sample, biological category, furthest sample groups based on the top PC dimensions retained, distance of each sample along a multidimensional line linking the most extreme samples, 0-1 scaled distance with direction based on expression level of most extreme groups.
#'   - information: data frame with category, number of variables (genes/proteins/others), method used ("PCA"), number of PCs kept, sum of variance kept, correlation between values extracted for samples in that biological category and their expression levels.
#'   - dimred_information: data frame with information needed for plot_dimensions(): biological category, number of PCs retained, percentage of variance described by PC1 and PC2, centroid, and slope of main axis of variance.
#'   When using method = "lda", the function returns a list with:
#'   - projection: data frame with sample, biological category, furthest sample groups based on the top PC dimensions retained, distance of each sample along a multidimensional line linking the most extreme samples, 0-1 scaled distance with direction based on expression level of most extreme groups.
#'   - information: data frame with category, number of variables (genes/proteins/others), method used ("PCA" or "LDA"), number of PCs and LDs kept, sum of variance kept, correlation between values extracted for samples in that biological category and their expression levels.
#'   - dimred_information: data frame with information needed for plot_dimensions(): biological category, number of PCs retained, percentage of variance described by PC1, PC2, LD1, and LD2, centroid, and slope of main axis of variance.
#' @export
dim_reduction <- function(data_input, method = c("scale","pca","lda"),threshold=0.5,lda_threshold=0.8, focus="group",lda_focus="group",correlation_method="spearman",remove_effect_from=NULL) {
  if (length(method) != 1) {
    stop("Please select exactly one method: 'scale', 'pca', 'lda'. Note that pca and lda will only run for categories with at least 5 variables/genes.")
  }
  expr = data_input$expr
  gene_meta = data_input$gene_meta
  sample_meta = data_input$sample_meta

  # Ensure rownames exist
  if (is.null(rownames(sample_meta))) {
    stop("rownames of the sample metadata must correspond to the column names from the expression matrix. \n")
  }

  # Reorder metadata to match columns of expr_mat
  if (!all(colnames(expr) %in% rownames(sample_meta))) {
    warning("Some samples from the expression matrix (column names) were not found in the sample metadata. \n")
  }

  sample_meta = sample_meta[which(sample_meta$sample%in%colnames(expr)), ]

  # Warn if the order was changed
  if (!all(rownames(sample_meta) == colnames(expr))) {
    warning("Sample metadata was reordered to match the column order from the expression matrix. \n")
  }

  sample_meta = sample_meta[match(colnames(expr),sample_meta$sample),]
  data_input2 = data_input

  if(length(remove_effect_from)>0){
      design = model.matrix(as.formula(paste("~ ",paste(remove_effect_from,collapse="+"))),data = sample_meta)
      fit = limma::lmFit(expr, design)
      expr = limma::residuals.MArrayLM(fit, expr)
      print(paste0("The expression data was adjusted based on the following: ",paste("~ ",paste(remove_effect_from,collapse="+"))))
      data_input2$expr = expr
  }
  method <- match.arg(method, choices = c("pca", "scale","lda"))


  if (method == "pca") {
    # run PCA method with threshold
  #  thres_input <- readline(prompt = "Enter variance threshold for PCA (e.g., 0.5 to select all dimensions which added together account to at least of 50% of the total variance): ")
  #  thres <- as.numeric(thres_input)
  #  if (is.na(thres) || thres <= 0 || thres > 1) {
  #    stop("Invalid threshold. Please enter a numeric value between 0 and 1.")
  #  }
    res <- pca_method(data_input2,threshold,focus,correlation_method)
  } else if (method == "scale") {
    # run normalized average expression method with min_avg_expr filter
    res <- scale_method(data_input2)
  } else if (method =="lda"){
    res <- lda_method(data_input2, threshold, lda_threshold,focus,lda_focus,correlation_method)

  }

  return(res)
}
