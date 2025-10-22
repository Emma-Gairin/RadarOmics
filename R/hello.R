# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

library(dplyr)
library(tidyr)
library(ggradar)
library(ggplot2)
library(patchwork)
library(MASS)
#resultnat <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist/vsd_ocellaris.csv",
 #                     "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist/fortest/sampleinfo_ocellaris.csv",
  #                    "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist/fortest/genelist_manuscript.csv")

#result2 <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/vsd_ocellaris_billy28.csv",
 #                     "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/sampleinfo_ocellaris_billy28.csv",
  #                    "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist/fortest/genelist_manuscript.csv")

#resultbilly <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/vsd_ocellaris_billy.csv",
 #                     "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/sampleinfo_ocellaris_billy.csv",
  #                    "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/genelist_manuscript.csv")


resultbilly <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/vsd_ocellaris_billy.csv",
                           "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/sampleinfo_ocellaris_billy.csv",
                           "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/genelist_manuscript.csv")


result <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/cyanea/vsd_juv_field.csv",
                           "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/cyanea/sample_information.csv",
                           "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/cyanea/gene_list.csv")


#resultbillynat <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/vsd_ocellaris_billy_natacha_batchcorrect.csv",
 #                          "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/sampleinfo_ocellaris_28_nat.csv",
  #                         "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/genelist_manuscript.csv")

#test=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("scale"))
method="lda"
threshold=0.35
lda_threshold  =0.9
lda_focus = "Natural"
focus = "substance_concentration"
test2=dim_reduction(
  result,
  method = method,
  threshold=threshold,lda_threshold=lda_threshold,focus=focus,lda_focus=lda_focus)

forpca = test2$pca$glycolysis
forpca$sample = rownames(forpca)
forpca = merge(forpca,result$sample_meta,by="sample")
ggplot(forpca,aes(x=PC1,y=PC2,color=group,size=Natural))+geom_point()


forpca = test2$lda$`HPI axis`
forpca$sample = rownames(forpca)
forpca = merge(forpca,result$sample_meta,by="sample")
ggplot(forpca,aes(x=LD1,y=LD2,color=group,size=Natural))+geom_point()



#test2=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("pca"))
#test2=dim_reduction(resultbillynat$expr,resultbillynat$gene_meta,resultbillynat$sample_meta, method = c("pca"))
table(test2$information$method)

category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
#category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
category_list
category_list_names = c(#"ossification",
 # "oxidative stress",                        "histone deacetylation","immunity",
  "olfactory",
                        #"stomach",#"digestion",                        "appetite_regulation",
                        "glycolysis","krebs_cycle",
                      #  "fatty_acid_metabolism",
                        "cholesterol_biosynthesis","TH_pathway_down",
                      "TH_pathway_up","corticoids","HPI_axis","neurotransmission")

category_list_names = c(#"ossification",
  "heat response","apoptosis",
                        "oxidative stress",    #                    "histone deacetylation",
                        "immunity",
                       # "olfactory",
                        #"vision",

                       # "stomach","digestion",                        "appetite_regulation",
                        "glycolysis","krebs_cycle","cholesterol_biosynthesis",
                        "fatty_acid_metabolism",
                        "TH_pathway_down",
                        "TH_pathway_up","corticoids","HPI_axis","neurotransmission",                        "telomere")

category_list_names = c("glycolysis","beta-oxidation","urea cycle","krebs",
                        "fatty acid synthesis","cholesterol biosynthesis","
                        phototransduction","HPT axis and thyroid gland","Thyroid hormone signalling",
                        "HPI axis","corticoids","hsp","immunity","ROS","phase1","phase2","phase3")

#category_list_names = c("ossification","olfactory","phototransduction","neurotransmission","appetite_regulation",
#                        "cholesterol_biosynthesis",
#  "glycolysis","krebs_cycle",
#    "fatty_acid_metabolism", "stomach","digestion",
# "TH_pathway",
# "corticoids")

#colour = read.csv("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_colour.csv")
category_list=cbind(category_list_names,c(1:length(category_list_names)))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
#radars =plot_radar(resultbilly,test2,category_list,colour)
radars =plot_radar(result,test2,category_list)

# include wrap plots in plot radar

# make option to export all the prelim plots also

#ordered_list <- c("s1","s2","s3","s4","s5","s6","s7")
#ordered_list <- c("s2","s3","s4","s5","s6","s7")

#ordered_list <- c("s2_low","s2_high","s3_low","s3_high","s4_low","s4_high","s5_low","s5_high",
#                  "s6_low","s6_high","s7_low","s7_high")

ordered_list = unique(result$sample_meta$group[order(result$sample_meta$Natural)])
#wrap_plots(radars[ordered_list], ncol = 2, nrow = 6)
wrap_plots(radars[ordered_list],ncol=5,nrow=4)

if(method=="lda"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_billy_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".pdf"),height=20,width=7)
}
if(method=="pca"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_billy_",method,"_",threshold,".pdf"),height=20,width=7)

}






result <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/stefano/vsd_stefano.csv",
                            "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/stefano/sampleinfo_radar.csv",
                            "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/stefano/genelist.csv")


#resultbillynat <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/vsd_ocellaris_billy_natacha_batchcorrect.csv",
#                          "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/sampleinfo_ocellaris_28_nat.csv",
#                         "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/genelist_manuscript.csv")

#test=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("scale"))
method="lda"
threshold=0.5
lda_threshold  =0.9
lda_focus = "group"
focus = "group"
test2=dim_reduction(
  result,
  method = method,
  threshold=threshold,lda_threshold=lda_threshold,focus=focus,lda_focus=lda_focus)


forpca = test2$pca$iridophore
forpca$sample = rownames(forpca)
forpca = merge(forpca,result$sample_meta,by="sample")
ggplot(forpca,aes(x=PC1,y=PC2,color=group))+geom_point()



forpca = test2$lda$iridophore
forpca$sample = rownames(forpca)
forpca = merge(forpca,result$sample_meta,by="sample")
ggplot(forpca,aes(x=LD1,y=LD2,color=group))+geom_point()



if(method=="pca"){
  test2=dim_reduction(
    result,
    method = method,
    threshold=threshold)
}


#test2=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("pca"))
#test2=dim_reduction(resultbillynat$expr,resultbillynat$gene_meta,resultbillynat$sample_meta, method = c("pca"))
table(test2$information$method)

category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
#category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
category_list
category_list_names = c(#"ossification",
  # "oxidative stress",                        "histone deacetylation","immunity",
  "olfactory",
  #"stomach",#"digestion",                        "appetite_regulation",
  "glycolysis","krebs_cycle",
  #  "fatty_acid_metabolism",
  "cholesterol_biosynthesis","TH_pathway_down",
  "TH_pathway_up","corticoids","HPI_axis","neurotransmission")

category_list_names = c(#"ossification",
  "heat response","apoptosis",
  "oxidative stress",    #                    "histone deacetylation",
  "immunity",
  # "olfactory",
  #"vision",

  # "stomach","digestion",                        "appetite_regulation",
  "glycolysis","krebs_cycle","cholesterol_biosynthesis",
  "fatty_acid_metabolism",
  "TH_pathway_down",
  "TH_pathway_up","corticoids","HPI_axis","neurotransmission",                        "telomere")

category_list_names = c("glycolysis","beta-oxidation","urea cycle","krebs",
                        "fatty acid synthesis","cholesterol biosynthesis","
                        phototransduction","HPT axis and thyroid gland","Thyroid hormone signalling",
                        "HPI axis","corticoids","hsp","immunity","ROS","phase1","phase2","phase3")

#category_list_names = c("ossification","olfactory","phototransduction","neurotransmission","appetite_regulation",
#                        "cholesterol_biosynthesis",
#  "glycolysis","krebs_cycle",
#    "fatty_acid_metabolism", "stomach","digestion",
# "TH_pathway",
# "corticoids")
category_list_names = c("visual","xanthophore","melanophore","iridophore","ossification","osmoregulation",

                        "thyroid","corticoids","appetite","glycolysis","PDH_complex",
                        "betaoxi","krebs","fatty","cholesterol",  "gastrointestinal","digestion")

#colour = read.csv("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_colour.csv")
category_list=cbind(category_list_names,c(1:length(category_list_names)))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
#radars =plot_radar(resultbilly,test2,category_list,colour)
radars =plot_radar(result,test2,category_list)

# include wrap plots in plot radar

# make option to export all the prelim plots also

#ordered_list <- c("s1","s2","s3","s4","s5","s6","s7")
#ordered_list <- c("s2","s3","s4","s5","s6","s7")

#ordered_list <- c("s2_low","s2_high","s3_low","s3_high","s4_low","s4_high","s5_low","s5_high",
#                  "s6_low","s6_high","s7_low","s7_high")

ordered_list = c("DMSO","CPF","MPI","T3IOP","CPF.T3IOP")
#wrap_plots(radars[ordered_list], ncol = 2, nrow = 6)
#wrap_plots(radars[ordered_list],ncol=5,nrow=4)
wrap_plots(radars[ordered_list],ncol=3,nrow=2)

if(method=="lda"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/stefano/radar_stefano_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".pdf"),height=20,width=7)
}
if(method=="pca"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/stefano/radar_stefano_",method,"_",threshold,".pdf"),height=20,width=7)

}






result <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/vsd_zebrafish.csv",
                      "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/sampleinfo.csv",
                      "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/genelist_zf_manuscript.csv")

#result$expr=result$expr[,-which(colnames(result$expr)%in%"12_DMSO_Kontrolle4_96h_Eppi-14_Index-D5")]

#result$sample_meta=result$sample_meta[-which(result$sample_meta$sample%in%"12_DMSO_Kontrolle4_96h_Eppi-14_Index-D5"),]


#result$expr=result$expr[,-which(colnames(result$expr)%in%"11_DMSO_Kontrolle3_96h_Eppi-13_Index-C5")]
#
#result$sample_meta=result$sample_meta[-which(result$sample_meta$sample%in%"11_DMSO_Kontrolle3_96h_Eppi-13_Index-C5"),]
#
#
#result$expr=result$expr[,-which(colnames(result$expr)%in%"5_Rotenone_LC10_96h_Eppi-6_Index-E4")]
#
#result$sample_meta=result$sample_meta[-which(result$sample_meta$sample%in%"5_Rotenone_LC10_96h_Eppi-6_Index-E4"),]


# remove 96h_96h because it's weird

remove=result$sample_meta$sample[which(result$sample_meta$time_hpe==96)]

result$expr=result$expr[,-which(colnames(result$expr)%in%remove)]#"5_Rotenone_LC10_96h_Eppi-6_Index-E4")]

result$sample_meta=result$sample_meta[-which(result$sample_meta$sample%in%remove),]#%"5_Rotenone_LC10_96h_Eppi-6_Index-E4"),]



#resultbillynat <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/vsd_ocellaris_billy_natacha_batchcorrect.csv",
#                          "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/sampleinfo_ocellaris_28_nat.csv",
#                         "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/genelist_manuscript.csv")

#test=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("scale"))


method="pca"
threshold=0.7
lda_threshold  =0.8
lda_focus = "substance_concentration"
focus = "time_hpe"
test2=dim_reduction(
  result,
  method = method,
  threshold=threshold,lda_threshold=lda_threshold,focus=focus,lda_focus=lda_focus)

#if(method=="pca"){
#  test2=dim_reduction(
#    result,
#    method = method,
#    threshold=threshold)
#}


#test2=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("pca"))
#test2=dim_reduction(resultbillynat$expr,resultbillynat$gene_meta,resultbillynat$sample_meta, method = c("pca"))
table(test2$information$method)

category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
#category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
category_list
category_list_names = c(#"ossification",
  # "oxidative stress",                        "histone deacetylation","immunity",
  "olfactory",
  #"stomach",#"digestion",                        "appetite_regulation",
  "glycolysis","krebs_cycle",
  #  "fatty_acid_metabolism",
  "cholesterol_biosynthesis","TH_pathway_down",
  "TH_pathway_up","corticoids","HPI_axis","neurotransmission")

category_list_names = c(#"ossification",
  "heat response","apoptosis",
  "oxidative stress",    #                    "histone deacetylation",
  "immunity",
  # "olfactory",
  #"vision",

  # "stomach","digestion",                        "appetite_regulation",
  "glycolysis","krebs_cycle","cholesterol_biosynthesis",
  "fatty_acid_metabolism",
  "TH_pathway_down",
  "TH_pathway_up","corticoids","HPI_axis","neurotransmission",                        "telomere")

category_list_names = c("glycolysis","beta-oxidation","urea cycle","krebs",
                        "fatty acid synthesis","cholesterol biosynthesis","
                        phototransduction","HPT axis and thyroid gland","Thyroid hormone signalling",
                        "HPI axis","corticoids","hsp","immunity","ROS","phase1","phase2","phase3")

#category_list_names = c("ossification","olfactory","phototransduction","neurotransmission","appetite_regulation",
#                        "cholesterol_biosynthesis",
#  "glycolysis","krebs_cycle",
#    "fatty_acid_metabolism", "stomach","digestion",
# "TH_pathway",
# "corticoids")
category_list_names = c("ROS",
                        "mitochondria",
                        "endoplasmic_reticulum",
                        "histone",
                        "angiogenesis",
                        "brain_development",
                        "ossification",
                        #"vasculogenesis",
                        "vision",
                        "pigmentation",

                      #  "glutathione",
                        #"osmoregulation",

                        "thyroid",
                      "corticoids",
                        #"appetite",
                        "glycolysis",
                        #"betaoxi",
                        "krebs",#"fatty",
                        "cholesterol")#,  "gastrointestinal")#,"digestion")

#colour = read.csv("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_colour.csv")
category_list=cbind(category_list_names,c(1:length(category_list_names)))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
#radars =plot_radar(resultbilly,test2,category_list,colour)
radars =plot_radar(result,test2,category_list)

# include wrap plots in plot radar

# make option to export all the prelim plots also

#ordered_list <- c("s1","s2","s3","s4","s5","s6","s7")
#ordered_list <- c("s2","s3","s4","s5","s6","s7")

#ordered_list <- c("s2_low","s2_high","s3_low","s3_high","s4_low","s4_high","s5_low","s5_high",
#                  "s6_low","s6_high","s7_low","s7_high")

unique(test2$projection$group)
ordered_list = c("DMSO_36_12","DMSO_48_24","DMSO_96_72",#"DMSO_96_96",
                 "Sorafenib_36_12","Sorafenib_48_24","Sorafenib_96_72",#"Sorafenib_96_96",
                 "Rotenone_36_12","Rotenone_48_24","Rotenone_96_72")
                 #,"Rotenone_96_96")
ordered_list=unique(test2$projection$group)[order(unique(test2$projection$group))]


#ordered_list = c(ordered_list[1:11],"DMSO_0_36_12",ordered_list[12:14],"DMSO_0_36_12",ordered_list[15:18])
#wrap_plots(radars[ordered_list], ncol = 2, nrow = 6)
#wrap_plots(radars[ordered_list],ncol=5,nrow=4)
wrap_plots(radars[ordered_list],ncol=5,nrow=3)

if(method=="lda"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/radar_zebrafish_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".pdf"),height=20,width=20)
  write.csv(test2$projection, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/projection_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".csv"))#,height=20,width=20)
  write.csv(test2$information, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/information_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".csv"))#,height=20,width=20)

}
if(method=="pca"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/radar_zebrafish_",method,"_",threshold,".pdf"),height=20,width=20)
  write.csv(test2$projection, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/projection_",method,"_",threshold,"_",focus,".csv"))#,height=20,width=20)
  write.csv(test2$information, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/zebrafish/information_",method,"_",threshold,"_",focus,".csv"))#,height=20,width=20)

}

# also add function to get correlation between categories
# add function to only display categories that show significant differences across groups
# and add function for index of dissimilarity between plots
# add function to perform radar plots on distinct groups distinctly like each dev stage independently of previous and post
# do a sensitivity analysis to lda threshold and normal threshold and plot correlation between shapes across both dimesnions, presenting shapes in corners.



# NATACHA

result <- import_data("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist_ncbi/vsd_ocellaris.csv",
                          "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist_ncbi/sampleinfo_ocellaris.csv",
                         "C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist_ncbi/genelist_manuscript.csv")

#test=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("scale"))
method="pca"
threshold=0.4
lda_threshold  =0.8
lda_focus = "substance"
focus = "group"
test2=dim_reduction(
  result,
  method = method,
  threshold=threshold,lda_threshold=lda_threshold,focus=focus,lda_focus=lda_focus)

plots=plot_dimensions(result,test2,colour="group")
test2$information[which(test2$information$category=="appetite"),]
head(test2$projection[which(test2$projection$category=="appetite"),]
)
plots$appetite


table(test2$information$method)

category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
category_list
category_list_names = c(
  "appetite","glycolysis","lactic","krebs","betaoxi",
                        "cholesterol","fatty","digestion","gastrointestinal",
                        "corticoids",
                        "thyroid","ossification","vision",
                        "pigmentation","melanophore","iridophore","xanthophore")
#colour = read.csv("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_colour.csv")
category_list=cbind(category_list_names,c(1:length(category_list_names)))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
#radars =plot_radar(resultbilly,test2,category_list,colour)
radars =plot_radar(result,test2,category_list)

# include wrap plots in plot radar

# make option to export all the prelim plots also

ordered_list <- c("s1","s2","s3","s4","s5","s6","s7")
#ordered_list <- c("s2","s3","s4","s5","s6","s7")

#ordered_list <- c("s2_low","s2_high","s3_low","s3_high","s4_low","s4_high","s5_low","s5_high",
#                  "s6_low","s6_high","s7_low","s7_high")

ordered_list = unique(result$sample_meta$group)
#wrap_plots(radars[ordered_list], ncol = 2, nrow = 6)
wrap_plots(radars[ordered_list],ncol=3,nrow=3)

if(method=="lda"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist_ncbi/radar_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".pdf"),height=7,width=7)
}
if(method=="pca"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_kallist_ncbi/radar_",method,"_",threshold,".pdf"),height=7,width=7)

}




Billy

# NATACHA

result <- resultbilly
#test=dim_reduction(result$expr,result$gene_meta,result$sample_meta, method = c("scale"))
method="pca"
threshold=0.75
lda_threshold  =0.8
lda_focus = "substance"
focus = "group"
test2=dim_reduction(
  result,
  method = method,
  threshold=threshold,lda_threshold=lda_threshold,focus=focus,lda_focus=lda_focus)

plots=plot_dimensions(result,test2,colour="group")
test2$information[which(test2$information$category=="appetite"),]
head(test2$projection[which(test2$projection$category=="appetite"),]
)
plots$appetite


table(test2$information$method)

category_list = cbind(unique(result$gene_meta$category),rev(c(1:length(unique(result$gene_meta$category)))))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
category_list
category_list_names = c(
  "appetite_regulation","glycolysis","krebs_cycle","beta_oxidation",
  "cholesterol_biosynthesis","fatty_acid_metabolism","digestion","stomach",
  "HPI_axis","corticoids","TH_pathway_up","TH_pathway_down",
  "ossification","vision", "heat response","protein","apoptosis","histone deacetylation")
#  "pigmentation","melanophore","iridophore","xanthophore")
#colour = read.csv("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_colour.csv")
category_list=cbind(category_list_names,c(1:length(category_list_names)))
colnames(category_list)=c("category","order")
category_list=as.data.frame(category_list)
#radars =plot_radar(resultbilly,test2,category_list,colour)
radars =plot_radar(result,test2,category_list)

# include wrap plots in plot radar

# make option to export all the prelim plots also

ordered_list <- c("s1","s2","s3","s4","s5","s6","s7")
#ordered_list <- c("s2","s3","s4","s5","s6","s7")

#ordered_list <- c("s2_low","s2_high","s3_low","s3_high","s4_low","s4_high","s5_low","s5_high",
#                  "s6_low","s6_high","s7_low","s7_high")

ordered_list = unique(result$sample_meta$group)
#wrap_plots(radars[ordered_list], ncol = 2, nrow = 6)
wrap_plots(radars[ordered_list],ncol=6,nrow=3)

if(method=="lda"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".pdf"),height=7,width=7)
  write.csv(test2$projection, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/projection_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".csv"))#,height=20,width=20)
  write.csv(test2$information, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/information_",method,"_",threshold,"_",lda_threshold,"_",focus,"_",lda_focus,".csv"))#,height=20,width=20)
}
if(method=="pca"){
  ggsave(paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/radar_",method,"_",threshold,".pdf"),height=7,width=7)
  write.csv(test2$projection, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/projection_",method,"_",threshold,"_",focus,".csv"))#,height=20,width=20)
  write.csv(test2$information, paste0("C:/Users/EMMA-GAIRIN/Documents/THESIS/Draft papers cyanea/radar method/ocellaris_dev_heat/information_",method,"_",threshold,"_",focus,".csv"))#,height=20,width=20)

}

