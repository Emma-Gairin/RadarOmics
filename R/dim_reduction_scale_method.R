#' Dimensional reduction based on 0-1 scaling.
#'
#' To be used as part of dim_reduction().
#'
#' For each gene category, computes average expression per sample,
#' normalizes to 0-1 range across samples,
#' calculates mean normalized expression per group,
#' identifies furthest groups by min and max mean expression,
#' returns per-sample values and metadata in PCA-like output format.
#'
#' @param data_input Expression data, sample information, and biological process list uploaded using import_data()
#' @return A list with:
#'   - projection: data frame with sample, normalised_distance (normalized avg expr), category, furthest_groups.
#'   - information: data frame with number of genes per category.
#' @export

scale_method <- function(data_input,focus="group") {
  # Align rows and columns
  expr = data_input$expr
  gene_meta = data_input$gene_meta
  sample_meta = data_input$sample_meta
  expr <- expr[intersect(gene_meta$gene, rownames(expr)), intersect(sample_meta$sample, colnames(expr)), drop = FALSE]

  catlist <- unique(gene_meta$category)
  if(length(which(colnames(sample_meta)%in%focus))==0){
    stop("focus (default: group) column not found in sample information table. Make sure to have a column named group or use argument focus.")
  }
  projection <- data.frame(sample = character(),
                                   normalised_distance = numeric(),
                                   category = character(),
                                   furthest_groups = character(),
                                   stringsAsFactors = FALSE)

  information <- data.frame(category = character(), n_variables = integer(), stringsAsFactors = FALSE)

  for (cat in catlist) {
    genes_in_cat <- gene_meta$gene[gene_meta$category == cat]
    genes_in_cat <- intersect(genes_in_cat, rownames(expr))
    n_variables <- length(genes_in_cat)
    information <- rbind(information, data.frame(category = cat, n_variables = n_variables, stringsAsFactors = FALSE))

    if (n_variables > 0) {
      sub_expr <- expr[genes_in_cat, , drop = FALSE]

      sub_expr_01 <- t(apply(sub_expr, 1, function(x) {
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      }))

      avg_expr <- colMeans(sub_expr_01, na.rm = TRUE)

      # Min-max normalization (0-1)
      norm_expr <- (avg_expr - min(avg_expr)) / (max(avg_expr) - min(avg_expr))

      df_samples <- data.frame(
        sample = names(norm_expr),
        normalised_distance = as.numeric(norm_expr),
        category = cat,
        stringsAsFactors = FALSE
      )

      # Merge with group info
      df_samples <- merge(df_samples, sample_meta[, c("sample", focus)], by = "sample")

      # Compute mean normalized expression per group
      group_means <- aggregate(normalised_distance ~ group, df_samples, mean)

      # Identify furthest groups by min and max mean normalized expression
      furthest_groups <- paste0(
        group_means$group[which.min(group_means$normalised_distance)],
        "-",
        group_means$group[which.max(group_means$normalised_distance)]
      )

      # Add furthest groups column for all samples in this category
      df_samples$furthest_groups <- furthest_groups

      # Combine into final output
      projection <- rbind(projection,
                                  df_samples[, c("sample",  "category", "furthest_groups","normalised_distance")])
    }
  }
  projection <- projection %>%
    dplyr::left_join(sample_meta[, c("sample", focus)], by = "sample")


  list(
    projection = projection,
    information = information
  )
}
