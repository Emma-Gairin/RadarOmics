#' Plot radars
#'
#' Plot radar charts following dim_reduction()
#'
#' @param data_input Expression data, sample information, and gene list uploaded using import_data()
#' @param dim_reduction_output Output from dim_reduction()
#' @param category_list Optional data frame with columns: 1) category names, 2) order on the radar plot. Default will follow the order of gene categories from input file used in import_data(). You do not have to include all categories originally provided in import_data(), and you can play around with the order in which categories are displayed to manually highlight relationships between biological processes.
#' @param colour Optional colour scheme choice for different groups of samples. For optimal visualisation quality, we recommend exporting the outputs as PDFs and manually editing the colour scheme and figure design with vector editing software.
#' @return A list of `ggradar` plots, one per group of samples.
#' @export
#'
#' @examples
#' plot_list=plot_radar(data_input,dim_reduction_output,category_list=my_category_list)

plot_radar=function(data_input,dim_reduction_output,radar_grouping=NULL,category_list=NA,colour_sample=NULL,colour_average=NULL,axis_label_size=1,radar_label_size=4,radar_label_position="center",width=2,height=2,control = NULL, scale_control_to_mean = 0.5) {
  # Check required packages
  if (!requireNamespace("ggradar",quietly=TRUE)) stop("Package 'ggradar' is required.",call.=FALSE)
  if (!requireNamespace("dplyr",quietly=TRUE)) stop("Package 'dplyr' is required.",call.=FALSE)
  if (!requireNamespace("tidyr",quietly=TRUE)) stop("Package 'tidyr' is required.",call.=FALSE)
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Package 'ggplot2' is required.",call.=FALSE)
  if(length(radar_grouping)<1){
    radar_grouping=colnames(dim_reduction_output$projection)[length(colnames(dim_reduction_output$projection))]
  }
  if(length(control)>0){
    i=0
    while(i<length(unique(dim_reduction_output$projection$category))){
      i=i+1
      catsel = unique(dim_reduction_output$projection$category)[i]
      dimsel = dim_reduction_output$projection[which(dim_reduction_output$projection$category%in%catsel),]
      mean_dmso = mean(dimsel$normalised_distance[which(dimsel[[radar_grouping]]==control)])
      rescaling = max(c(1-mean_dmso),mean_dmso)/(abs(scale_control_to_mean-mean_dmso))
      dimsel$normalised_distance = dimsel$normalised_distance/rescaling
      mean_dmso = mean(dimsel$normalised_distance[which(dimsel[[radar_grouping]]==control)])

      dimsel$normalised_distance = dimsel$normalised_distance + abs(scale_control_to_mean-mean_dmso)-scale_control_to_mean
      mult = scale_control_to_mean/(max(c(abs(max(dimsel$normalised_distance))),abs(c(min(dimsel$normalised_distance)))))
      dimsel$normalised_distance = dimsel$normalised_distance*mult + scale_control_to_mean

      dim_reduction_output$projection[which(dim_reduction_output$projection$category%in%catsel),] = dimsel
    }

  }

  if(length(radar_grouping)<1){
    radar_grouping=colnames(dim_reduction_output$projection)[length(colnames(dim_reduction_output$projection))]
  }
  projection=dim_reduction_output$projection
  projection = merge(projection,
                     data_input$sample_meta[ , !(names(data_input$sample_meta) %in% intersect(names(projection), names(data_input$sample_meta)[-which(names(data_input$sample_meta)=="sample")]))],
                     by="sample")

  plot_list=list()
  list_target=unique(projection[[radar_grouping]])

  if(radar_label_position == "top"){radar_label_position = 1.8}
  if(radar_label_position == "center"){radar_label_position = 0}

  # If category_list not provided,default to order in data_input$gene_meta
  if (length(category_list)==1) {
    categories=unique(data_input$gene_meta$category)
    category_list=data.frame(
      category=categories,
      order=seq_along(categories)
    )
  } else {
    # ensure order column is numeric
    category_list$order=as.numeric(category_list$order)
  }
  for (target in list_target) {
    radar_target=projection %>%
      dplyr::filter(!!sym(radar_grouping) == !!target) %>%
      dplyr::select(sample,category,normalised_distance)

    reshaped_radar=radar_target %>%
      tidyr::pivot_wider(names_from=category,values_from=normalised_distance)

    nsample_per_group=nrow(reshaped_radar)
    reshaped_radar2=rbind(reshaped_radar,c(sample = target,colMeans(reshaped_radar[,-1])))

    meanradar=rbind(rep(0,ncol(reshaped_radar2)),
                      rep(1,ncol(reshaped_radar2)),
                      as.data.frame(reshaped_radar2))

    # Convert to numeric
    meanradar[-1]=lapply(meanradar[-1],function(x) as.numeric(as.character(x)))

    # Set sample factor levels
    if(as.character(target)%in%c("0","1")){print("Error. Target / group names should not be 0 or 1")}
    desired_order=c("0","1",reshaped_radar$sample,target)
    meanradar$sample=factor(meanradar$sample,levels=desired_order)

    # Custom colors

    if(length(colour_sample)==0){
      custom_colors=c("transparent","transparent",rep("#d6d6d6",nsample_per_group),"black")
      if(length(colour_average)>0){
        custom_colors=c("transparent","transparent",rep("#d6d6d6",nsample_per_group),colour_average$colour[which(colour_average[[radar_grouping]]==target)])
      }

    }
    if(length(colour_sample)>0){
      meanradar_colour = merge(meanradar,colour_sample,by="sample")
      custom_colors=c("transparent","transparent",meanradar_colour$colour,"black")
      if(length(colour_average)>0){
        custom_colors=c("transparent","transparent",meanradar_colour$colour,colour_average$colour[which(colour_average[[radar_grouping]]==target)])
      }
    }

    names(custom_colors)=desired_order

    # Determine column ordering based on category_list
    categories_meanradar=category_list$category[category_list$category %in% colnames(meanradar)]
    colnames_ordering=categories_meanradar[order(category_list$order[category_list$category %in% colnames(meanradar)])]

    meanradar_organised=meanradar[,c("sample",colnames_ordering)]

    # Generate radar plot
    radar_plot=suppressMessages(ggradar(meanradar_organised,
                         group.colours=custom_colors,
                         group.line.width=1,
                         group.point.size=0,
                         grid.line.width=0.5,
                         grid.mid=0,
                         axis.label.size=axis_label_size,
                         grid.label.size=0,
                         legend.position="NONE") +
      xlim(-width, width) +   # increase upper limit to create space
      ylim(min(-height,-radar_label_position), max(height,radar_label_position)) +   # increase upper limit to create space
      annotate("text", x = 0, y = radar_label_position, label = target,
               size = radar_label_size, fontface = "bold"))
      #annotate("text",x=0,y=1.5,label=group,size=radar_label_size,fontface="bold")

    plot_list[[target]]=radar_plot
  }

  return(plot_list)
}
