#' Print PCA/LDA plots and if <3 dimensions are retained, display main axis onto which each sample is projected.
#'
#' @param data_input Expression data, sample information, and gene list uploaded using import_data()
#' @param dim_reduction_output Output from dim_reduction()
#' @param colour Define colour of plotted data points based on column name from the sample information matrix (default = "group")
#' @param shape Define shape of plotted data points based on column name from the sample information matrix (default = NA)
#' @return A list of plots displaying the coordinates of each sample for each biological category

plot_dimensions=function(data_input,dim_reduction_output,colour="group",shape=NA) {
  # Check required packages
  if (!requireNamespace("dplyr",quietly=TRUE)) stop("Package 'dplyr' is required.",call.=FALSE)
  if (!requireNamespace("tidyr",quietly=TRUE)) stop("Package 'tidyr' is required.",call.=FALSE)
  if (!requireNamespace("ggplot2",quietly=TRUE)) stop("Package 'ggplot2' is required.",call.=FALSE)
  plot_list=list()

  for(cat in unique(dim_reduction_output$information$category)){
    reduction_method = dim_reduction_output$information$method[which(dim_reduction_output$information$category==cat)][1]
    if(reduction_method=="PCA"){
      reduction_dim = dim_reduction_output$information$num_pcs[which(dim_reduction_output$information$category%in%cat)][1]
      pc1 = dim_reduction_output$information$pc1[which(dim_reduction_output$information$category%in%cat)][1]
      pc2 = dim_reduction_output$information$pc2[which(dim_reduction_output$information$category%in%cat)][1]

      forpca = dim_reduction_output$pca[[cat]]
      forpca$sample = rownames(forpca)
      forpca = merge(forpca,data_input$sample_meta,by="sample")

      forpca[[colour]] <- as.factor(forpca[[colour]])
      if(length(shape)>1){      forpca[[shape]]  <- as.factor(forpca[[shape]])}

      centroid = dim_reduction_output$information$centroid[which(dim_reduction_output$information$category%in%cat)]
      maxvariancedirection = dim_reduction_output$information$maxvariancedirection[which(dim_reduction_output$information$category%in%cat)]
      if(length(centroid)==0){
        if(length(shape)>1){
          plot = ggplot(forpca,aes(x=PC1,y=PC2,color=!!sym(colour),shape=!!sym(shape)))+geom_point()+
            xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
            ggtitle(paste0(cat," | PC retained:",reduction_dim))+coord_equal()
        }
        if(length(shape)==1){
          plot = ggplot(forpca,aes(x=PC1,y=PC2,color=!!sym(colour)))+geom_point()+
            xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
            ggtitle(paste0(cat," | PC retained:",reduction_dim))+coord_equal()
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
        if(length(shape)>1){
        plot = ggplot(forpca,aes(x=PC1,y=PC2,color=!!sym(colour),shape=!!sym(shape)))+geom_point()+
          geom_point() +
            geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
          xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
          ggtitle(paste0(cat," | PC retained:",reduction_dim))+coord_equal()
        }
        if(length(shape)==1){
          plot = ggplot(forpca,aes(x=PC1,y=PC2,color=!!sym(colour)))+geom_point()+
            geom_point() +
            geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
            xlab(paste0("PC1 (",round(100*pc1),"%)"))+ylab(paste0("PC2 (",round(100*pc2),"%)"))+
            ggtitle(paste0(cat," | PC retained:",reduction_dim))+coord_equal()
        }
        plot_list[[cat]]=plot

      }
    }

      if(reduction_method=="LDA"){
        reduction_dim = dim_reduction_output$information$num_lds[which(dim_reduction_output$information$category%in%cat)][1]
        forlda = dim_reduction_output$lda[[cat]]
        forlda$sample = rownames(forlda)
        forlda = merge(forlda,data_input$sample_meta,by="sample")
        forlda[[colour]] <- as.factor(forlda[[colour]])
        if(length(shape)>1){      forlda[[shape]]  <- as.factor(forlda[[shape]])}

        centroid = dim_reduction_output$information$centroid[which(dim_reduction_output$information$category%in%cat)]
        maxvariancedirection = dim_reduction_output$information$maxvariancedirection[which(dim_reduction_output$information$category%in%cat)]
        if(lda_focus!=focus){
          if(length(shape)>1){
            plot = ggplot(forlda,aes(x=LD1,y=LD2,color=!!sym(colour),shape=!!sym(shape)))+geom_point()+
              geom_point() +
              ggtitle(paste0(cat))+
              ggtitle(paste0(cat," | LD retained:",reduction_dim))
          }
          if(length(shape)==1){

          plot = ggplot(forlda,aes(x=LD1,y=LD2,color=!!sym(colour)))+geom_point()+
            geom_point() +
            ggtitle(paste0(cat))+
            ggtitle(paste0(cat," | LD retained:",reduction_dim))
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
          if(length(shape)>1){

          plot = ggplot(forlda,aes(x=LD1,y=LD2,color=!!sym(colour),shape=!!sym(shape)))+geom_point()+
            geom_point() +
            geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
            ggtitle(paste0(cat))+coord_equal()+
          ggtitle(paste0(cat," | LD retained:",reduction_dim))
          }
          if(length(shape)==1){
            plot = ggplot(forlda,aes(x=LD1,y=LD2,color=!!sym(colour)))+geom_point()+
              geom_point() +
              geom_abline(slope = slope, intercept = intercept, color = "black", lwd = 1)+
              ggtitle(paste0(cat))+coord_equal()+
              ggtitle(paste0(cat," | LD retained:",reduction_dim))
          }
          plot_list[[cat]]=plot

        }
    }
  }
  return(plot_list)
}
