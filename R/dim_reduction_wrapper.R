#' Dimensional reduction of expression counts using 0-1 scaled expression levels, PCA, or LDA
#'
#' For each category of variables (e.g., group of genes related to a given biological process), reduces the dimensionality across variables to a single value for each sample.
#'
#' If using method="scale", the function returns the average value of the 0-1 scaled variable within each category for each sample
#'
#' If using method="pca", the function uses a Principal Component Analysis (PCA) on the expression matrix to derive the coordinates of the samples for each group. The top PC dimensions are retained based on a variance threshold that can be modified in the function (default: threshold = 0.5). The axis which maximises the variance in sample-averaged PCA coordinates for the user-defined focus (default: focus = "group") is defined. Each sample is projected onto this axis. The most extreme samples are set to have values of 0 and 1, with lower projection values corresponding to the side of the axis showing lower expression values.
#'
#' If using method="lda", the function uses a Principal Component Analysis (PCA) on the expression matrix to derive the coordinates of the samples for each group. The top PC dimensions are retained based on a variance threshold that can be modified in the function (default: threshold = 0.5). A Linear Discriminant Analysis (LDA) is used to increase the separation between samplles based on a user-defined lda-focus (default: lda_focus = "group"). The top LD dimensions are retained based on a variance threshold that can be modified in the function (default: lda_threshold = 0.8). The axis which maximises the variance in sample-averaged LDA coordinates for the user-defined focus (default: focus = "group") is defined. Each sample is projected onto this axis. The most extreme samples are set to have values of 0 and 1, with lower projection values corresponding to the side of the axis showing lower expression values.
#'
#' @param data_input Expression data, sample information, and biological process list uploaded using import_data()
#' @return
#' @export
dim_reduction <- function(data_input, method = c("scale","pca","lda"),threshold=0.5,lda_threshold=0.8, focus="group",lda_focus="group") {
  if (length(method) != 1) {
    stop("Please select exactly one method: 'scale', 'pca', 'lda'. Note that pca and lda will only run for categories with at least 5 variables/genes.")
  }
  expr = data_input$expr
  gene_meta = data_input$gene_meta
  sample_meta = data_input$sample_meta

  method <- match.arg(method, choices = c("pca", "scale","lda"))


  if (method == "pca") {
    # run PCA method with threshold
  #  thres_input <- readline(prompt = "Enter variance threshold for PCA (e.g., 0.5 to select all dimensions which added together account to at least of 50% of the total variance): ")
  #  thres <- as.numeric(thres_input)
  #  if (is.na(thres) || thres <= 0 || thres > 1) {
  #    stop("Invalid threshold. Please enter a numeric value between 0 and 1.")
  #  }
    res <- pca_method(data_input, threshold,focus)
  } else if (method == "scale") {
    # run normalized average expression method with min_avg_expr filter
    res <- scale_method(data_input)
  } else if (method =="lda"){
    res <- lda_method(data_input, threshold, lda_threshold,focus,lda_focus)

  }

  return(res)
}
