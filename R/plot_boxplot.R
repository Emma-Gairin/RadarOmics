#' Plot boxplots for data inspection
#'
#' Produce boxplots of the values extracted for each sample across categories, grouping samples based on default = "group" or user defined criteria.
#'
#' @param data_input Counts data, sample information, and biological feature list uploaded using import_data()
#' @param dim_reduction_output Output from dim_reduction()
#' @param stat_choice Choice of "nonparametric" (Kruskal + Dunn), or "parametric" (ANOVA + Tukey's HSD test). Default = "nonparametric".
#' @param pvalue_adjustment Choice of p-value adjustment methods for non-parametric testing. Same options as for dunn.test() from library(dunn.test). Default = "bh" (Benjamini-Hochberg)
#' @param pvalue_threshold Choice of p-value significance threshold (default = 0.05)
#' @param boxplot_grouping Choice of sample grouping for boxplot visualisation and statistical testing. Default = "group" from data_input$sample_meta.
#' @param boxplot_order Manual modification of the order of the groups of boxplots. Must be provided as a vector, e.g., boxplot_order = c("s1","s2",...).
#' @param colour_palette Define colour palette (default: a mix of palettes from the "wesanderson" package).
#' @return A list of plots displaying the values yielded by dim_reduction() for each sample in each biological category. One plot per category.

plot_boxplot=function(data_input,dim_reduction_output,stat_choice = "nonparametric",pvalue_adjustment="bh",pvalue_threshold=0.05,boxplot_focus="group",boxplot_order = NULL,colour_palette = c(wes_palette("AsteroidCity1"),
                                                                                                           wes_palette("Chevalier1"),wes_palette("Darjeeling2"),
                                                                                                                 wes_palette("Darjeeling1"),
                                                                                                                 wes_palette("GrandBudapest1"),wes_palette("Moonrise3"),wes_palette("Moonrise2"))) {
  # Check required packages
  if (!requireNamespace("dplyr",quietly=TRUE)) stop("Package 'dplyr' is required.",call.=FALSE)
  if (!requireNamespace("tidyr",quietly=TRUE)) stop("Package 'tidyr' is required.",call.=FALSE)
  if (!requireNamespace("ggplot2::ggplot2",quietly=TRUE)) stop("Package 'ggplot2::ggplot2' is required.",call.=FALSE)
  if (!requireNamespace("wesanderson",quietly=TRUE)) stop("Package 'wesanderson' is required.",call.=FALSE)
  plot_list=list()

  for(cat in unique(dim_reduction_output$projection$category)){
    stat_choice2 = 1
    projection_cat=dim_reduction_output$projection[which(dim_reduction_output$projection$category==cat),]
    if(length(intersect(boxplot_focus,colnames(projection_cat)))==0){
      boxplot_focus = colnames(projection_cat)[ncol(projection_cat)]
    }
    projection_cat = merge(projection_cat,
                       data_input$sample_meta[ , !(names(data_input$sample_meta) %in% intersect(names(projection_cat), names(data_input$sample_meta)[-which(names(data_input$sample_meta)=="sample")]))],
                       by="sample")
    if(length(boxplot_order)>0){
      projection_cat[,which(colnames(projection_cat)%in%boxplot_focus)] = factor(projection_cat[,which(colnames(projection_cat)%in%boxplot_focus)],levels = boxplot_order)
    }

    if(stat_choice == "parametric"){
      set.seed(123)
      if (!requireNamespace("multcompView",quietly=TRUE)) stop("Package 'multcompView' is required.",call.=FALSE)

      projection_cat$test = projection_cat[,which(colnames(projection_cat)%in%boxplot_focus)]
      aov_test = aov(normalised_distance  ~ test, data = projection_cat)
      tukey_test = TukeyHSD(aov_test)
      projection_cat_letters = multcompView::multcompLetters4(aov_test, tukey_test)
      projection_cat_letters = as.data.frame.list(projection_cat_letters$test)
      projection_cat_letters$test = rownames(projection_cat_letters)
      projection_cat = merge(projection_cat,projection_cat_letters[,c("test","Letters")],by="test")
      if(length(unique(projection_cat$Letters))==1){
        boxplot = ggplot2::ggplot(projection_cat,aes(x=!!rlang::sym(boxplot_focus),y=normalised_distance,fill = !!rlang::sym(boxplot_focus)))+geom_boxplot()+
          theme_bw()+xlab(paste0(boxplot_focus))+ylab("dim_reduction() value")+theme(legend.position = "none")+
          scale_fill_manual(values=colour_palette)+labs(title=paste0(cat),subtitle=c(" No significant difference"))
      }
      if(length(unique(projection_cat$Letters))>1){

      boxplot = ggplot2::ggplot(projection_cat,aes(x=!!rlang::sym(boxplot_focus),y=normalised_distance,fill = !!rlang::sym(boxplot_focus)))+geom_boxplot()+
        theme_bw()+xlab(paste0(boxplot_focus))+ylab("dim_reduction() value")+
        scale_fill_manual(values=colour_palette)+ggtitle(paste0(cat))+theme(legend.position = "none")+
        geom_text(aes(x = !!rlang::sym(boxplot_focus), y = 1.1, label = Letters), size = 3.5)
      }
      }
    if(stat_choice == "nonparametric"){
      if (!requireNamespace("FSA",quietly=TRUE)) stop("Package 'FSA' is required.",call.=FALSE)

      projection_cat$test = projection_cat[,which(colnames(projection_cat)%in%boxplot_focus)]
      kruskal_test = kruskal.test(x=projection_cat$normalised_distance,g = projection_cat$test)
      if(kruskal_test$p.value<pvalue_threshold){
        dunn_test = dunn.test(x=projection_cat$normalised_distance, g=projection_cat$test, method = pvalue_adjustment)
        dunn_test_result = cbind(dunn_test$comparisons,dunn_test$P.adjusted)
        dunn_test_result = as.data.frame(dunn_test_result)
        significant = setNames(
          dunn_test_result$V2 < pvalue_threshold,      # TRUE if significant
          gsub(" ", "", dunn_test_result$V1) # remove spaces
        )

        projection_cat_letters = multcompView::multcompLetters(significant)

        projection_cat_letters = as.data.frame(projection_cat_letters$Letters)
        projection_cat_letters$test = rownames(projection_cat_letters)
        colnames(projection_cat_letters)=c("Letters","test")
        projection_cat = merge(projection_cat,projection_cat_letters[,c("Letters","test")],by="test")

        boxplot = ggplot2::ggplot(projection_cat,aes(x=!!rlang::sym(boxplot_focus),y=normalised_distance,fill = !!rlang::sym(boxplot_focus)))+geom_boxplot()+
          theme_bw()+xlab(paste0(boxplot_focus))+ylab("dim_reduction() value")+
          scale_fill_manual(values=colour_palette)+ggtitle(paste0(cat))+theme(legend.position = "none")+
          geom_text(aes(x = !!rlang::sym(boxplot_focus), y = 1.1, label = Letters), size = 3.5)


      }
      if(kruskal_test$p.value>pvalue_threshold){
        stat_choice2 = NULL
      }

    }
    if(min(length(stat_choice),length(stat_choice2))==0){

      boxplot = ggplot2::ggplot(projection_cat,aes(x=!!rlang::sym(boxplot_focus),y=normalised_distance,fill = !!rlang::sym(boxplot_focus)))+geom_boxplot()+
        theme_bw()+xlab(paste0(boxplot_focus))+ylab("dim_reduction() value")+theme(legend.position = "none")+
        scale_fill_manual(values=colour_palette)+ggtitle(paste0(cat))
      if(length(stat_choice2)==0){
        boxplot = ggplot2::ggplot(projection_cat,aes(x=!!rlang::sym(boxplot_focus),y=normalised_distance,fill = !!rlang::sym(boxplot_focus)))+geom_boxplot()+
          theme_bw()+xlab(paste0(boxplot_focus))+ylab("dim_reduction() value")+theme(legend.position = "none")+
          scale_fill_manual(values=colour_palette)+labs(title=paste0(cat),subtitle=c(" No significant difference"))
      }
    }
      plot_list[[cat]]=boxplot
    }

    return(plot_list)
  }

