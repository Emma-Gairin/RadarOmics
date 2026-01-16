#' PCA/LDA plots
#'
#' Produce PCA/LDA plots for PC1/2 or LD1/2. If 2+ dimensions are required to meet the variance threshold used in dim_reduction(), display main axis onto which each sample is projected.
#'
#' @param data_input Counts data, sample information, and molecular feature list uploaded using import_data()
#' @param dim_reduction_output Output from dim_reduction()
#' @param colour Define colour of plotted data points based on column name from the sample information matrix (default = "group")
#' @param shape Define shape of plotted data points based on column name from the sample information matrix (default = NA)
#' @param point_size Define size of plotted data points (default = 3)
#' @param colour_palette Define colour palette (default: a mix of palettes from the "wesanderson" package)
#' @return A list of plots displaying the PC1/PC2 or LD1/LD2 coordinates of each sample for each biological category and main axis of variance when 2 dimensions are retained by dim_reduction(). One plot per category.


plot_dimensions=function(data_input,dim_reduction_output,colour="group",shape="",point_size=3,colour_palette = c(wes_palette("AsteroidCity1"),
                                                                                                    wes_palette("Chevalier1"),wes_palette("Darjeeling2"),
                                                                                                    wes_palette("Darjeeling1"),
                                                                                                    wes_palette("GrandBudapest1"),wes_palette("Moonrise3"),wes_palette("Moonrise2"))) {
  # Check required packages
  if (!requireNamespace("dplyr",quietly=TRUE)) stop("Package 'dplyr' is required.",call.=FALSE)
  if (!requireNamespace("tidyr",quietly=TRUE)) stop("Package 'tidyr' is required.",call.=FALSE)
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Package 'ggplot2' is required.",call.=FALSE)
  if (!requireNamespace("wesanderson",quietly=TRUE)) stop("Package 'wesanderson' is required.",call.=FALSE)
  plot_list=list()

  for(cat in unique(dim_reduction_output$dimred_information$category)){
    reduction_method = dim_reduction_output$dimred_information$method[which(dim_reduction_output$dimred_information$category==cat)][1]
    if(reduction_method=="PCA"){
      reduction_dim = dim_reduction_output$dimred_information$num_pcs[which(dim_reduction_output$dimred_information$category%in%cat)][1]
      pc1 = dim_reduction_output$dimred_information$pc1[which(dim_reduction_output$dimred_information$category%in%cat)][1]
      pc2 = dim_reduction_output$dimred_information$pc2[which(dim_reduction_output$dimred_information$category%in%cat)][1]

      forpca = dim_reduction_output$pca[[cat]]
      forpca$sample = rownames(forpca)
      forpca = merge(forpca,data_input$sample_meta,by="sample")

      forpca[[colour]] <- as.factor(forpca[[colour]])
      if(nchar(shape)>0){      forpca[[shape]]  <- as.factor(forpca[[shape]])}

      centroid = dim_reduction_output$dimred_information$centroid[which(dim_reduction_output$dimred_information$category%in%cat)]
      maxvariancedirection = dim_reduction_output$dimred_information$maxvariancedirection[which(dim_reduction_output$dimred_information$category%in%cat)]
      if(length(centroid)==0){
        if(nchar(shape)>0){
          plot = ggplot(forpca,aes(x=PC1,y=PC2,color=factor(!!rlang::sym(colour)),shape=factor(!!rlang::sym(shape))))+geom_point(size=point_size)+
            xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
            ggtitle(paste0(cat," | PC retained: ",reduction_dim))+coord_equal()+theme_bw()+scale_colour_manual(values = colour_palette)
        }
        if(nchar(shape)==0){
          plot = ggplot(forpca,aes(x=PC1,y=PC2,color=factor(!!rlang::sym(colour))))+geom_point(size=point_size)+
            xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
            ggtitle(paste0(cat," | PC retained: ",reduction_dim))+coord_equal()+theme_bw()+scale_colour_manual(values = colour_palette)
        }


        plot_list[[cat]]=plot

      }
      if(length(centroid)>0){
        centroid_PC1=centroid[1]
        centroid_PC2=centroid[2]
        slope_PC1=maxvariancedirection[1]
        slope_PC2=maxvariancedirection[2]
        slope=slope_PC2/slope_PC1
        intercept=centroid_PC1-slope*centroid_PC1
        centroid_df = data.frame(PC1 = centroid_PC1, PC2 = centroid_PC2)
        if(nchar(shape)>0){
        plot = ggplot(forpca,aes(x=PC1,y=PC2,color=factor(!!rlang::sym(colour)),shape=factor(!!rlang::sym(shape))))+geom_point(size=point_size)+
          geom_point(size=point_size) +
            geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
          xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
          ggtitle(paste0(cat," | PC retained: ",reduction_dim))+coord_equal()+theme_bw()+scale_colour_manual(values = colour_palette)
        }
        if(nchar(shape)==0){
          plot = ggplot(forpca,aes(x=PC1,y=PC2,color=factor(!!rlang::sym(colour))))+geom_point(size=point_size)+
            geom_point(size=point_size) +
            geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
            xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
            ggtitle(paste0(cat," | PC retained: ",reduction_dim))+coord_equal()+theme_bw()+scale_colour_manual(values = colour_palette)
        }
        plot_list[[cat]]=plot

      }
    }

      if(reduction_method=="LDA"){
        reduction_dim = dim_reduction_output$dimred_information$num_lds[which(dim_reduction_output$dimred_information$category%in%cat)][1]
        ld1 = dim_reduction_output$dimred_information$ld1[which(dim_reduction_output$dimred_information$category%in%cat)][1]
        ld2 = dim_reduction_output$dimred_information$ld2[which(dim_reduction_output$dimred_information$category%in%cat)][1]
        forlda = dim_reduction_output$lda[[cat]]
        forlda$sample = rownames(forlda)
        forlda = merge(forlda,data_input$sample_meta,by="sample")
        forlda[[colour]] <- as.factor(forlda[[colour]])
        if(nchar(shape)>0){      forlda[[shape]]  <- as.factor(forlda[[shape]])}

        centroid = dim_reduction_output$dimred_information$centroid[which(dim_reduction_output$dimred_information$category%in%cat)]
        maxvariancedirection = dim_reduction_output$dimred_information$maxvariancedirection[which(dim_reduction_output$dimred_information$category%in%cat)]
          if(nchar(shape)>0){
            plot = ggplot2::ggplot(forlda,ggplot2::aes(x=LD1,y=LD2,color=factor(!!rlang::sym(colour)),shape=factor(!!rlang::sym(shape))))+ggplot2::geom_point(size=point_size)+
              geom_point(size=point_size) +
              xlab(paste0("LD1 (",round(100*ld1),"%)"))+ylab(paste0("LD2 (",round(100*ld2),"%)"))+
              ggtitle(paste0(cat," | LD retained: ",reduction_dim))

          if(nchar(shape)==0){

          plot = ggplot(forlda,aes(x=LD1,y=LD2,color=factor(!!rlang::sym(colour))))+geom_point(size=point_size)+
            geom_point(size=point_size) +
            xlab(paste0("LD1 (",round(100*ld1),"%)"))+ylab(paste0("LD2 (",round(100*ld2),"%)"))+
            ggtitle(paste0(cat," | LD retained: ",reduction_dim))+theme_bw()+scale_colour_manual(values = colour_palette)
          }
        }
        if(length(centroid)>0){
          centroid_LD1=centroid[1]
          centroid_LD2=centroid[2]
          slope_LD1=maxvariancedirection[1]
          slope_LD2=maxvariancedirection[2]
          slope=slope_LD2/slope_LD1
          intercept=centroid_LD1-slope*centroid_LD1
          centroid_df = data.frame(LD1 = centroid_LD1, LD2 = centroid_LD2)
          if(nchar(shape)>0){

          plot = ggplot(forlda,aes(x=LD1,y=LD2,color=factor(!!rlang::sym(colour)),shape=factor(!!rlang::sym(shape))))+geom_point(size=point_size)+
            geom_point(size=point_size) +
            geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
            xlab(paste0("LD1 (",round(100*ld1),"%)"))+ylab(paste0("LD2 (",round(100*ld2),"%)"))+
            ggtitle(paste0(cat))+coord_equal()+theme_bw()+scale_colour_manual(values = colour_palette)+
          ggtitle(paste0(cat," | LD retained: ",reduction_dim))
          }
          if(nchar(shape)==0){
            plot = ggplot2::ggplot(forlda,aes(x=LD1,y=LD2,color=factor(!!rlang::sym(colour))))+geom_point(size=point_size)+
geom_point(size=point_size) +
              geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
              ggtitle(paste0(cat))+coord_equal()+theme_bw()+scale_colour_manual(values = colour_palette)+
              xlab(paste0("LD1 (",round(100*ld1),"%)"))+ylab(paste0("LD2 (",round(100*ld2),"%)"))+
              ggtitle(paste0(cat," | LD retained: ",reduction_dim))
          }
          plot_list[[cat]]=plot

        }
    }
  }
  return(plot_list)
}
