#' Dimensional reduction based on 0-1 scaling.
#'
#' To be used as part of dim_reduction().
#'
#' For each biological category, computes average expression per sample,
#' normalizes to 0-1 range across samples,
#' calculates mean normalized counts per group,
#' identifies furthest groups by min and max mean counts,
#' returns per-sample values and metadata in PCA-like output format.
#'
#' @param data_input Counts data, sample information, and biological process list uploaded using import_data()
#' @return A list with:
#'   - projection: data frame with sample, normalised_distance, category, furthest_groups.
#'   - information: data frame with number of features per category.
#' @export

scale_method <- function(data_input,focus="group") {
  # Align rows and columns
  counts = data_input$counts
  feature_meta = data_input$feature_meta
  sample_meta = data_input$sample_meta
  counts <- counts[intersect(feature_meta$feature, rownames(counts)), intersect(sample_meta$sample, colnames(counts)), drop = FALSE]

  catlist <- unique(feature_meta$category)
  if(length(which(colnames(sample_meta)%in%focus))==0){
    stop("focus (default: group) column not found in sample information table. Make sure to have a column named group or use argument focus.")
  }
  projection <- data.frame(sample = character(),
                                   normalised_distance = numeric(),
                                   category = character(),
                                   furthest_groups = character(),
                                   stringsAsFactors = FALSE)

  information <- data.frame(category = character(), n_features = integer(), stringsAsFactors = FALSE)

  for (cat in catlist) {
    features_in_cat <- feature_meta$feature[feature_meta$category == cat]
    features_in_cat <- intersect(features_in_cat, rownames(counts))
    n_features <- length(features_in_cat)
    information <- rbind(information, data.frame(category = cat, n_features = n_features, stringsAsFactors = FALSE))

    if (n_features > 0) {
      sub_counts <- counts[features_in_cat, , drop = FALSE]

      sub_counts_01 <- t(apply(sub_counts, 1, function(x) {
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      }))

      avg_counts <- colMeans(sub_counts_01, na.rm = TRUE)

      # Min-max normalization (0-1)
      norm_counts <- (avg_counts - min(avg_counts)) / (max(avg_counts) - min(avg_counts))

      df_samples <- data.frame(
        sample = names(norm_counts),
        normalised_distance = as.numeric(norm_counts),
        category = cat,
        stringsAsFactors = FALSE
      )

      # Merge with group info
      df_samples <- merge(df_samples, sample_meta[, c("sample", focus)], by = "sample")

      # Compute mean normalized counts per group
      group_means <- aggregate(normalised_distance ~ group, df_samples, mean)

      # Identify furthest groups by min and max mean normalized counts
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
  projection <-
    dplyr::left_join(projection,sample_meta[, c("sample", focus)], by = "sample")


  list(
    projection = projection,
    information = information
  )
}
