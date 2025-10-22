#' Plot radar charts from PCA dim_reduction_output
#'
#' @param data_input Expression data, sample information, and gene list uploaded using import_data()
#' @param dim_reduction_output Output from dim_reduction()
#' @param category_list Optional data frame with columns: 1) category names,2) order on the radar plot. Default will follow the order of gene categories from input file used in import_data()
#'
#' @return A list of `ggradar` plots,one per group
#' @export
#'
#' @examples
#' plot_list=plot_radar(data_input,dim_reduction_output,category_list=my_category_list)
plot_radar_groups=function(data_input,dim_reduction_output,category_list=NULL) {
  # Check required packages
  if (!requireNamespace("ggradar",quietly=TRUE)) stop("Package 'ggradar' is required.",call.=FALSE)
  if (!requireNamespace("dplyr",quietly=TRUE)) stop("Package 'dplyr' is required.",call.=FALSE)
  if (!requireNamespace("tidyr",quietly=TRUE)) stop("Package 'tidyr' is required.",call.=FALSE)
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Package 'ggplot2' is required.",call.=FALSE)

  projection=dim_reduction_output$projection
  plot_list=list()
  list_groups=unique(projection$group)

  # If category_list not provided,default to order in data_input$gene_meta
  if (is.null(category_list)) {
    categories=unique(data_input$gene_meta$category)
    category_list=data.frame(
      category=categories,
      order=seq_along(categories)
    )
  } else {
    # ensure order column is numeric
    category_list$order=as.numeric(category_list$order)
  }

  for (group in list_groups) {
    radar_group=projection %>%
      dplyr::filter(group == !!group) %>%
      dplyr::select(sample,category,normalised_distance)

    reshaped_radar=radar_group %>%
      tidyr::pivot_wider(names_from=category,values_from=normalised_distance)

    nsample_per_group=nrow(reshaped_radar)
    reshaped_radar2=rbind(reshaped_radar,c(group,colMeans(reshaped_radar[,-1])))

    meanradar=rbind(rep(0,ncol(reshaped_radar2)),
                    rep(1,ncol(reshaped_radar2)),
                    as.data.frame(reshaped_radar2))

    # Convert to numeric
    meanradar[-1]=lapply(meanradar[-1],function(x) as.numeric(as.character(x)))

    # Set sample factor levels
    desired_order=c("0","1",reshaped_radar$sample,group)
    meanradar$sample=factor(meanradar$sample,levels=desired_order)

    # Custom colors
    custom_colors=c("transparent","transparent",rep("#9f9f9fff",nsample_per_group),"black")
    names(custom_colors)=desired_order

    # Determine column ordering based on category_list
    categories_meanradar=category_list$category[category_list$category %in% colnames(meanradar)]
    colnames_ordering=categories_meanradar[order(category_list$order[category_list$category %in% colnames(meanradar)])]

    meanradar_organised=meanradar[,c("sample",colnames_ordering)]

    # Generate radar plot
    radar_plot=ggradar(meanradar_organised,
                       group.colours=custom_colors,
                       group.line.width=1,
                       group.point.size=0,
                       grid.line.width=0.5,
                       grid.mid=0,
                       axis.label.size=2.5,
                       grid.label.size=0,
                       legend.position="NONE") +
      annotate("text",x=0,y=0,label=group,size=6,fontface="bold")

    plot_list[[group]]=radar_plot
  }

  return(plot_list)
}
