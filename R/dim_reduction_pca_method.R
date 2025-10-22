#' For each gene category,calculation of distance between samples based on top dimensions from Principal Component Analyses
#'
#' Perform PCA on the expression level of genes in a given gene category,
#' identify the number of PCs corresponding to a certain threshold of variance,
#' identify the axis of maximum variance between the average coordinates of each group,
#' project the samples onto this axis in the multidimensional space,
#' determine which side of the axis is associated with lower values based on average gene expression in 2 most extreme groups.
#'
#' @param expr Gene expression matrix (genes as rownames,samples as column names,cell fill as gene expression level). This gene expression matrix should be normalised (e.g.,variance-stabilised).
#' @param sample_meta Sample metadata (one row per sample. column names of this table should include "sample","group". within the sample column,the names should match the column names of the expression matrix)
#' @param gene_meta Gene metadata (one row per gene. column names of this table should include "gene","category". within the gene column,the names should feature within the row names of the expression matrix)
#' @param threshold threshold of cumulative variance used to identify the number of PCs to keep.
#' @return A list with:
#'   - projection: data frame with distance of each sample along a multidimensional line linking the most extreme samples,category,and furthest group pair.
#'   - information: data frame with category,number of PCs kept,sum of variance kept,and number of genes considered. Distances obtained for the PCA
#' @export
pca_method = function(data_input,threshold,focus){
  expr = data_input$expr
  gene_meta = data_input$gene_meta
  sample_meta = data_input$sample_meta
  gene_list_prep = intersect(rownames(expr),unique(unique(gene_meta$gene)))
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
    num_pcs = integer(),
    sum_variance_kept = numeric(),
    n_genes = integer(),
    expr_pca_correlation = numeric(),
    flipped = logical(),
    stringsAsFactors = FALSE
  )
  pca = list()
  # for each gene category,
  for (cat in catlist){
    # extract the rows with genes belonging to the category
    genes_in_cat = gene_meta$gene[gene_meta$category == cat]
    filtered_expr = expr[rownames(expr) %in% genes_in_cat,,drop = FALSE]
    # calculate how many genes we retrieve
    n_genes = nrow(filtered_expr)
    # for very low number of genes,we recommend not to run the PCA and to use scaled expression instead
    if (n_genes > 4){
      # PCA on transposed matrix with samples as rows and genes as columns
      filtered_expr = filtered_expr[apply(filtered_expr, 1, sd) != 0,]

      sample_pca = prcomp(t(filtered_expr),scale. = TRUE)

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
      num_pcs = min(which(cumulative_variance > threshold))
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



      flipped = FALSE
      if (avg_diff < 0){
        coordinate_mainaxis$distance2 = -coordinate_mainaxis$distance

        # Normalize distance to 0-1 after flipping if needed
        min_val = min(coordinate_mainaxis$distance2)
        max_val = max(coordinate_mainaxis$distance2)

        coordinate_mainaxis$normalised_distance = (coordinate_mainaxis$distance2 - min_val)/ (max_val - min_val)

        flipped = TRUE
      }

      if (avg_diff >= 0){
        coordinate_mainaxis$distance2 = coordinate_mainaxis$distance

        # Normalize distance to 0-1 after flipping if needed
        min_val = min(coordinate_mainaxis$distance2)
        max_val = max(coordinate_mainaxis$distance2)

        coordinate_mainaxis$normalised_distance = (coordinate_mainaxis$distance2 - min_val)/ (max_val - min_val)

      }
      # Correlation between normalized coordinate and average gene expression per sample in category
      avg_expr_per_sample = colMeans(filtered_expr)
      avg_expr_per_sample = avg_expr_per_sample[coordinate_mainaxis$sample]
      correlation = suppressWarnings(cor.test(coordinate_mainaxis$normalised_distance,avg_expr_per_sample,use = "complete.obs",method="spearman"))
      corr_val = correlation$estimate
      p_value = correlation$p.value

      # Save summary info including correlation and flipped status
      information = rbind(
        information,
        data.frame(
          category = cat,
          method="PCA",
          num_pcs = num_pcs,
          centroid = centroid,
          maxvariancedirection = maxvariancedirection,
          sum_variance_kept_pcs = sum_var_kept,
          pc1 = pc1,
          pc2 = pc2,
          n_genes = n_genes,
          expr_pca_correlation_spearman_rho = corr_val,
          expr_pca_correlation_pvalue = p_value,
          flipped = flipped,
          stringsAsFactors = FALSE
        )
      )
      rownames(information) = NULL
      projection = rbind(projection,
                         data.frame(sample = coordinate_mainaxis$sample,
                                    category = cat,
                                    furthest_groups = furthest_groups,
                                    distance = coordinate_mainaxis$distance,
                                    normalised_distance = coordinate_mainaxis$normalised_distance,
                                    stringsAsFactors = FALSE))

    }
  }
  projection <- projection %>%
    dplyr::left_join(sample_meta[, c("sample", "group")], by = "sample")

  list(projection = projection,
       information = information,
       pca=pca)
}
