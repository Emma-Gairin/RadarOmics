#' Dimensional reduction based on PCAs.
#'
#' To be used as part of dim_reduction().
#' For each biological category, extraction of a single value per sample based on top dimensions from Principal Component Analyses
#'
#' Perform PCA on the counts matrix of genes/proteins/others from a biological category.
#' identify the number of PCs corresponding to a user-defined threshold of variance,
#' if PC1 meets the threshold, extract the PC1 coordinate of each sample.
#' If not, identify the main axis of variance between the average PC coordinates of each group (or other user-defined sample grouping scheme),
#' project the samples onto this axis in the multidimensional space,
#' determine which side of the axis is associated with lower values based on average gene expression in 2 most extreme groups,
#' scale sample values from 0 to 1, setting 0 to be the on the side of lower values.
#'
#' @param data_input Expression data, sample information, and biological process list uploaded using import_data()
#' @param threshold threshold of cumulative variance used to identify the number of PCs to keep,
#' @param focus grouping of samples used to determine the main axis of variance (default = "group")
#' @param correlation_method statistical test used to determine linear correlation between values extracted for each sample and mean scaled gene expression in the category (default = "spearman").
#' @return A list with:
#'   - projection: data frame with sample, biological category, furthest sample groups based on the top PC dimensions retained, distance of each sample along a multidimensional line linking the most extreme samples, 0-1 scaled distance with direction based on expression level of most extreme groups.
#'   - information: data frame with category, number of variables (genes/proteins/others), method used ("PCA"), number of PCs kept, sum of variance kept, correlation between values extracted for samples in that biological category and their expression levels.
#'   - dimred_information: data frame with information needed for plot_dimensions(): biological category, number of PCs retained, percentage of variance described by PC1 and PC2, centroid, and slope of main axis of variance.
#' @export
pca_method = function(data_input,pca_threshold,focus,correlation_method,pca_scale){
  expr = data_input$expr
  gene_meta = data_input$gene_meta
  sample_meta = data_input$sample_meta
  gene_list_prep = intersect(rownames(expr),unique(gene_meta$gene))
  expr = expr[gene_list_prep,sample_meta$sample,drop = FALSE]

  catlist = unique(gene_meta$category)

  # preparation of the output of the function
  projection = data.frame(sample = character(),
                          distance = numeric(),
                          category = character(),
                          furthest_groups = character(),
                          stringsAsFactors = FALSE)

  # some extra information about the run -- number of Principal Components considered based on user choice of variance threshold
  information = data.frame(
    category = character(),
    method= character(),
    num_pcs = integer(),
    sum_variance_kept = numeric(),
    n_variables = integer(),
    expr_pca_correlation = numeric(),
    flipped = logical(),
    stringsAsFactors = FALSE
  )
  pca = list()
  dimred_information = data.frame(    category = character(),
                                   method= character(),
                                   num_pcs = integer(),
                                   pc1 = integer(),
                               pc2 = integer(),
                               centroid = integer(),
                               maxvariancedirection = integer())
  # for each gene category,
  for (cat in catlist){
    # extract the rows with genes belonging to the category
    genes_in_cat = gene_meta$gene[gene_meta$category == cat]
    filtered_expr = expr[rownames(expr) %in% genes_in_cat,,drop = FALSE]
    # calculate how many genes we retrieve
    n_variables = nrow(filtered_expr)
    # for very low number of genes,we recommend not to run the PCA and to use scaled expression instead
    if (n_variables > 4){
      # PCA on transposed matrix with samples as rows and genes as columns
      filtered_expr = filtered_expr[apply(filtered_expr, 1, sd) != 0,]

      sample_pca = prcomp(t(filtered_expr),scale. = pca_scale)

      prep_mean_group = as.data.frame(sample_pca$x)
      pca[[cat]] = prep_mean_group
      prep_mean_group$sample = rownames(prep_mean_group)

      # Merge with sample_meta to match samples with groups
      prep_mean_group = merge(prep_mean_group,sample_meta[,c("sample",focus)],by = "sample")
      rownames(prep_mean_group) = prep_mean_group$sample
      prep_mean_group$sample = NULL

      # Calculate the variance and cumulative variance corresponding to each dimension of the PCA
      explained_variance = sample_pca$sdev^2
      proportion_variance = explained_variance/sum(explained_variance)
      cumulative_variance = cumsum(proportion_variance)

      # Number of PCs to keep
      num_pcs = min(which(cumulative_variance > pca_threshold))
      sum_var_kept = sum(proportion_variance[1:num_pcs])
      pc1 = proportion_variance[1]
      pc2 = proportion_variance[2]
      prep_mean_group$focus = prep_mean_group[,focus]
      # Mean coordinate per group for all PCs kept
      mean_group = aggregate(.~focus,prep_mean_group[,c(paste0("PC",1:num_pcs),"focus")],mean)

      # Function to project sample coordinates onto main axis of variance
      projection_multidim = function(samplecoord,centroid,maxvariancedirection){
        samplecentroid = as.numeric(samplecoord-centroid)
        projectionmainaxis = sum(samplecentroid*maxvariancedirection)/sum(maxvariancedirection^2)
        projectedvalue = centroid+projectionmainaxis*maxvariancedirection
        projectedvalue
      }

      if (num_pcs==1){
        coordinate_mainaxis = data.frame(distance = prep_mean_group[,"PC1"])
        rownames(coordinate_mainaxis) = rownames(prep_mean_group)
        centroid=NA
        maxvariancedirection=NA
      } else{
        centroid = colMeans(mean_group[,paste0("PC",1:num_pcs)])
        covariancematrix = cov(mean_group[,paste0("PC",1:num_pcs)])
        eigen_data_input = eigen(covariancematrix)
        maxvariancedirection = eigen_data_input$vectors[,1]
        projected_points = t(apply(prep_mean_group[,paste0("PC",1:num_pcs)],1,projection_multidim,
                                   centroid = centroid,maxvariancedirection = maxvariancedirection))
        coordinate_mainaxis = as.data.frame(projected_points%*%maxvariancedirection)
        colnames(coordinate_mainaxis) = "distance"
      }

      coordinate_mainaxis$sample = rownames(coordinate_mainaxis)
      coordinate_mainaxis = merge(coordinate_mainaxis,sample_meta[,c("sample",focus)],by = "sample")
      coordinate_mainaxis$focus = coordinate_mainaxis[,focus]
      mean_group_dim = aggregate(. ~ focus,coordinate_mainaxis[,c("distance","focus")],mean)
      furthest_groups = paste0(
        mean_group_dim$focus[which.min(mean_group_dim$distance)],
        "-",
        mean_group_dim$focus[which.max(mean_group_dim$distance)]
      )

      # Get average expression in the two furthest groups for genes in category
      group1 = mean_group_dim$focus[which.min(mean_group_dim$distance)]
      group2 = mean_group_dim$focus[which.max(mean_group_dim$distance)]

      expr_group1 = filtered_expr[,sample_meta$sample[sample_meta[,focus] == group1],drop = FALSE]
      expr_group2 = filtered_expr[,sample_meta$sample[sample_meta[,focus] == group2],drop = FALSE]

      avg_expr_group1 = rowMeans(expr_group1)
      avg_expr_group2 = rowMeans(expr_group2)
      avg_diff = mean(avg_expr_group2 - avg_expr_group1)


      flipped=FALSE

      if (avg_diff < 0){
        coordinate_mainaxis$distance2 = -coordinate_mainaxis$distance
        flipped=TRUE

        # Normalize distance to 0-1 after flipping if needed
        min_val = min(coordinate_mainaxis$distance2)
        max_val = max(coordinate_mainaxis$distance2)

        coordinate_mainaxis$normalised_distance = (coordinate_mainaxis$distance2 - min_val)/ (max_val - min_val)

      }

      if (avg_diff >= 0){
        coordinate_mainaxis$distance2 = coordinate_mainaxis$distance
        flipped=FALSE

        # Normalize distance to 0-1 after flipping if needed
        min_val = min(coordinate_mainaxis$distance2)
        max_val = max(coordinate_mainaxis$distance2)

        coordinate_mainaxis$normalised_distance = (coordinate_mainaxis$distance2 - min_val)/ (max_val - min_val)

      }
      # Correlation between normalized coordinate and average gene expression per sample in category

      sub_expr_01 = t(apply(filtered_expr, 1, function(x) {
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      }))

      avg_expr = colMeans(sub_expr_01, na.rm = TRUE)

      # Min-max normalization (0-1)
      norm_expr = (avg_expr - min(avg_expr)) / (max(avg_expr) - min(avg_expr))

      avg_expr_per_sample = norm_expr
      avg_expr_per_sample = avg_expr_per_sample[coordinate_mainaxis$sample]
      correlation = suppressWarnings(cor.test(coordinate_mainaxis$normalised_distance,avg_expr_per_sample,use = "complete.obs",method=correlation_method))
      corr_val = correlation$estimate
      p_value = correlation$p.value

      # Save summary information
      information = rbind(
        information,
        data.frame(
          category = cat,
          n_variables = n_variables,
          method="PCA",
          num_pcs = num_pcs,
          sum_variance_kept_pcs = sum_var_kept,
          expr_pca_correlation = corr_val,
          expr_pca_correlation_pvalue = p_value,
          flipped=flipped,
          stringsAsFactors = FALSE
        )
      )

      rownames(information) = NULL

      dimred_information = rbind(dimred_information,
                              data.frame(category = cat,
                                         method="PCA",
                                         num_pcs = num_pcs,
                                         pc1 = pc1,
                                         pc2 = pc2,
                                         centroid = centroid,
                                         maxvariancedirection = maxvariancedirection))
      rownames(dimred_information) = NULL

      projection = rbind(projection,
                         data.frame(sample = coordinate_mainaxis$sample,
                                    category = cat,
                                    furthest_groups = furthest_groups,
                                    distance = coordinate_mainaxis$distance,
                                    normalised_distance = coordinate_mainaxis$normalised_distance,
                                    stringsAsFactors = FALSE))

    }
  }
  projection = projection %>%
    dplyr::left_join(sample_meta[, unique(c("sample", focus))], by = "sample")

  list(projection = projection,
       information = information,
       dimred_information = dimred_information,
       pca=pca)
}
