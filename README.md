# RadarOmics

RadarOmics is an R package for summarising **'omics datasets** through **dimensional reduction** and visualising the output across samples, treatments, and biological processes with **radar plots**.
It is designed to handle **complex experimental designs** with >2 treatments, and can handle nested designs (_e.g._, multiple developmental stages and chemical exposure treatments). It works with gene expression, protein expression, or similar tabular datasets.

We provide examples of outputs with simple datasets in **Implementation** below.

 We developed this method as part of a _landscape transcriptomics_ study with a complex design: juvenile and adult fish samples from two developmental stages and 18 field sites ([Gairin et al. (bioRxiv)](https://doi.org/10.1101/2025.08.05.668635)).

 ![Radar plots from the landscape transcriptomics study](images/example.png)
*Figure 1: Example output of the package - following manual editing on vector graphics software*

#
## Detailed description
RadarOmics summarises the expression profiles of genes or other molecules within predefined biological categories using dimensional reduction analysis and radar plot visualisations.

Prior to using the package, three data frames must be prepared in CSV (Comma delimited) format: 1) a count matrix normalised according to the preferred dimensional reduction method, e.g., variance-stabilising transformation of gene counts prior to running Principal Component Analysis; 2) sample information including the name of the samples and their corresponding group; 3) biological information with the gene/protein/etc. identifiers and their corresponding user-defined categories based on *e.g.,* manually curated lists, Gene Ontology terms, KEGG pathways. See examples of dataset structure [below](#implementation).

---
The three main functions of the package are:

**import_data()** - used to upload three data frames with **read.csv()** count matrix, sample information, and biological information.

**dim_reduction()** - performs dimensional reduction based on the subsetted counts matrix for each biological category (_i.e,._ keeping only rows corresponding to genes/proteins/others of interest) and yields a single value per sample and per biological category. There are three options of dimensional reduction:
- method = **"scale"**,
- method = **"pca"**,
- method = **"lda"**

**plot_radar()**  - used to generate radar plots displaying the values from **dim_reduction** for each sample and each gene category. One radar plot is produced for each group of samples (default grouping is the "group" column of the sample information data frame, unless otherwise specified using _radar_grouping_). The order of categories in the radar plot follows the order of gene categories in the biological information file unless a different order is specified by the user.

Users can mix-and-match methods and manually create the table to feed into **plot_radar()**, in particular if some biological categories have too few samples to be analysed with method = **"pca"** or **"lda"**. Similarly, users can also run **dim_reduction()** on multiple types of datasets (_e.g.,_ combining RNAseq counts with metabolomics or with phenotype information), or add phenotypic information to the **dim_reduction()** output, before visualising the results in the radarplot with **plot_radar()**.

While we provide a few options for text size and colour coding as part of **plot**, we recommend exporting the output of **plot_radar()** as a PDF with **ggsave("radar.pdf")** and manually editing the radars using vector graphics software to obtain publication-quality visuals (we use [Inkscape](https://inkscape.org/)). For instance, colour shading can be added as a background to highlight biological categories belonging to similar processes (_e.g.,_ highlighting energy metabolism, endocrine processes, _etc._).


---
### How does dim_reduction() work?

- method = **"scale"** lineary scales the expression level of each gene/protein/other for each biological category to a value between 0 and 1 (0: sample with lowest expression level, 1: highest). The average scaled expression level across all genes/proteins/others for each biological category is then calculated for each sample. The value obtained for each sample is then linearly scaled again from 0 to 1 with the sample showing lowest average value assigned 0 and highest value 1. This yields one value between 0 and 1 for each sample and biological category.

Methods **"pca"** and **"lda"** operate similarly, and a sketch of the main steps is provided below ([Figure 2](#visual-package-description)).
- method = **"pca"** starts with generating a **Principal Component Analysis** (**prcomp(,scale=TRUE)** from base R package _stats_) based on the expression level of each gene/protein/other for each biological category. The minimum number of PC dimensions  accounting for a user-defined percentage of the total variance is selected (___threshold___ = *e.g.,* 0.25, 0.50, 0.75). If PC1 represents a higher portion of the variance than __threshold__ when inspecting a given biological category, the PC1 coordinate of each sample will be used. If more than 1 dimension is necessary to reach __threshold__, the average coordinates of all samples from each **group** for these top PC dimensions are calculated (thus, the choice of how to **group** samples together has an importance on the final result. By default, the package uses the column "group" but this can be modified using the argument ___focus___ in **dim_reduction()). The first eigenvector of the covariance matrix obtained from these averages defines the main axis of variance across groups. This main axis is drawn through the multidimensional space for all top PC dimensions selected. Each sample is projected onto this main axis in multidimensional space, yielding points along a segment bounded by the two most extreme samples. The two extrememost samples are assigned values of 0 and 1. To determine which sample has a value of 0 or 1, the expression levels of the two sites with the most extreme average projections along the segment are compared, and the extrememost sample closest to the site with lowest average expression level is set to have a value of 0, while the opposite extremity is set to 1.
  
- method = **"lda"** starts with generating a Principal Component Analysis (**pca = prcomp(,scale=TRUE)** from base R package _stats_) based on the expression level of each gene/protein/other for each biological category. The minimum number of PC dimensions (**num_pcs**) accounting for a user-defined percentage of the total variance is selected (___threshold___ = *e.g.,* 0.25, 0.50, 0.75). If PC1 represents a higher portion of the variance than __threshold__ when inspecting a given biological category, the PC1 coordinate of each sample will be used. If more than 1 dimension is necessary to reach __threshold__, the coordinates of each sample across these top dimensions are used in a Linear Discriminant Analysis (LDA) with the function **MASS::lda(pca$x[,1:num_pcs], grouping = lda_focus**. This LDA maximises the variance between samples based on a user-defined ___lda_focus___ in **dim_reduction()**, which by default is "group" but can be set to a different column from the sample information table. Once the LDA is performed, the minimum number of LD dimensions accounting for a user-defined percentage of the total variance is selected (___lda_threshold___ = *e.g.,* 0.8, 0.9). The average coordinates of all samples from each **group** for these top LD dimensions are calculated (thus, the choice of how to group samples together has an importance on the final result. By default, the package uses the column "group" but this can be modified using the argument ___focus___ in **dim_reduction()). The first eigenvector of the covariance matrix obtained from these averages defines the main axis of variance across groups. This main axis is drawn through the multidimensional space for all top LD dimensions selected. Each sample is projected onto this main axis in multidimensional space, yielding points along a segment bounded by the two most extreme samples. The two extrememost samples are assigned values of 0 and 1. To determine which sample has a value of 0 or 1, the expression levels of the two sites with the most extreme average projections along the segment are compared, and the extrememost sample closest to the site with lowest average expression level is set to have a value of 0, while the opposite extremity is set to 1. 

***
### Visual package description
For a given biological category, after running a PCA or LDA and selecting the top dimensions, one value is retained per sample. If a single dimension is necessary to meet a user-defined threshold of variance, the value is derived from the coordinate of each sample along the first dimension.

If more than 1 dimension is necessary, the samples are projected along the main axis of variance across the groups of samples. The method is illustrated here with 2 PC dimensions:
 ![Extracting a value for each sample using PCA](images/method.png)
*Figure 2: Methodological sketch describing steps from __dim_reduction()__, here using method = "pca".*

***
### Recommendations
- method = **"scale"** is recommended only in the case of biological categories with few genes/proteins/others (< 5 -10) when expression levels are consistent or positively correlated across the category
- method = **"pca"** is recommended when using the package with simple experimental designs (_e.g.,_ one dimensional developmental series with a few stages, see [example #1](#example-1---method--pca-and-method--scale)) or with complex experimental designs where users want to preserve the variance across all samples (particularly during data exploration or when the effect of various co-acting variables on the samples is similar in extent or depends on the biological category of interest).
- method = **"lda"** is recommended when using the package with complex experimental designs where one variable acting on the samples has a strong influence that is not the main focus of the study (_e.g.,_ developmental timeseries with a few stages and multiple treatments: in that case, on a PCA, the signal of the developmental stage obscures that of the treatment and so using _lda_focus = "treatment"_ would better disentangle the effect of the treatment on the samples, see [example #2](#example-2---method--lda)).

### Notes
RadarOmics is intended to facilitate data visualisation across many biological categories and samples, and while offering multiple analytical options, it is best used in combination with other approaches to validate the results.
In particular, for methods "pca" and "lda" we recommend testing multiple _threshold_ and _lda_threshold_ values, inspecting the output of **dim_reduction()**, producing multiple options of radar plots, and cross-checking results for each biological category with other visualisations (_e.g.,_ heatmaps for each biological category) before making a final choice.

We provide various data inspection solutions when using method = **"pca"** or **"lda"**,
- __$information__, output from __dim_reduction()__, provides the number of PC and LD dimensions retained for each category
- __$information__, output from __dim_reduction()__, provides the linear correlation between the projected coordinates from PCAs and LDAs and the average scaled gene expression of each sample for each category (in detail - 0-1 scaling of the expression across samples, average of all genes/proteins/others in biological category for each sample, 0-1 scaling of this value based on most extreme samples). The correlation is is calculated with a user-selected statistical test (_correlation_method_ argument in __dim_reduction()__), such as Pearson or Spearman (default: "spearman"),
- __$dimred_information__, outputs from __dim_reduction()__, provide the coordinates of the samples for each PCA and LDA generated. PC1/2 or LD1/2, along with the main axis of variance used to derive a value for each sample (only if the number of dimensions retained is > 1), can be plotted using **plot_dimension()**.



---

## Installation

The development version can be directly installed from GitHub using **remotes::install_github()**.

```r
# Install remotes if not already installed
install.packages("remotes")

# Install RadarOmics from GitHub
remotes::install_github("Emma-Gairin/RadarOmics", auth_token = "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")
```

---
## Implementation

## Example #1 - method = "pca" and method = "scale" 
Here we go through an example pipeline using method = **"pca"** and method = **"scale"** from **dim_reduction()** to summarise the gene expression profile, for a pre-defined set of biological categories, of samples from different groups.

We use the RNAseq data from the 7-stage developmental series of the false clownfish _Amphiprion ocellaris_ (from [Roux et al. (2023)](https://doi.org/10.1016/j.celrep.2023.112661)). There are 21 samples in total (3 samples per stage x 7 stages).

### Load the package
```r
# load the package
library(RadarOmics)
```

### Import the data
We are supplying:
- a variance-stabilisation transformed count table (vsd) obtained using DESEq2 with genes as rows, samples as columns
- sample information with two columns: samples and their grouping (here, developmental stage, from stage 1 to stage 7)
- gene list with two columns: genes and their categories (here, various biological categories, _e.g.,_ glycolysis, Krebs cycle, phototransduction). This gene list is derived from a manually curated database ([Herrera et al., 2025](https://doi.org/10.1002/jez.b.23299)].

```r
# import data
data_input = import_data(expr_path = "vsd_ocellaris.csv", sample_meta_path = "sampleinfo_ocellaris.csv", gene_meta_path = "genelist_ocellaris.csv")
```
#
Expression data (or other tabular data), normalised for PCA use. For gene expression data, we recommend VSD normalisation with DESEq2.
```r
head(data_input$expr[,1:10])
```
|               | SRR7610156| SRR7610157| SRR7610162| SRR7610144| SRR7610145| SRR7610163| SRR7610146| SRR7610147| SRR7610148| SRR7610149|
|:--------------|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
|YP_001054867.1 |   13.99266|   14.36672|   14.34720|   14.11891|   14.44253|   14.08794|   14.27059|   14.18767|   14.59782|   15.71519|
|YP_001054868.1 |   13.57260|   13.89591|   13.76110|   13.66079|   13.99851|   13.51833|   13.67022|   13.61725|   14.24787|   15.25302|
|YP_001054869.1 |   18.23788|   18.97548|   18.62184|   18.31808|   18.79498|   18.51865|   18.40833|   18.70388|   18.56178|   19.28712|
|YP_001054870.1 |   15.80561|   16.85047|   16.49806|   16.33397|   16.62497|   16.48044|   16.17988|   16.61900|   16.64477|   17.49150|
|YP_001054871.1 |   10.60595|   11.21210|   11.32445|   11.69303|   11.67610|   11.26584|   11.17176|   11.72576|   11.77252|   12.31938|
|YP_001054872.1 |   15.80708|   16.16411|   15.73523|   16.09713|   16.08820|   15.68903|   16.00612|   16.11608|   16.22958|   17.62626|

*Note that samples are columns, genes are rows.*
#
Sample information with columns "sample" and "group". One radar plot per "group" will be generated.
```r
head(data_input$sample_meta)
```
|sample     |group |
|:----------|:-----|
|SRR7610156 |s1    |
|SRR7610157 |s1    |
|SRR7610162 |s1    |
|SRR7610144 |s2    |
|SRR7610145 |s2    |
|SRR7610163 |s2    |

*Note that the package will look for a column named "group" unless otherwise specified as arguments in **dim_reduction()**.*
#
Gene information with columns "gene" and "category". Users can manually reshuffle and filter categories before plotting the output of the package using the radar plot.
```r
data_input$gene_meta
```
|gene           |category |
|:--------------|:--------|
|XP_023142913.1 |appetite |
|XP_023142914.1 |appetite |
|XP_023120868.1 |appetite |
...
|XP_054861428.1 |vision   |
|XP_054861429.1 |vision   |
|XP_023135802.2 |vision   |

*Note that the package needs a column named "category".*
#
### Dimensional reduction - PCA
Once the dataset is uploaded, we can run the PCA and extract reduced coordinates from each sample and each biological category based on top PC dimensions representing e.g., 50 % of variance (default _threshold = 0.5_).
```r
dim_reduction_output = dim_reduction(
  data_input,
  method = "pca")
```
#
**dim_reduction()** yields multiple objects.
By using scaling, PCA, or LDA, it obtains a value between 0 and 1 for each sample.
#
$projection provides, for each sample and biological category, the two furthest groups (e.g., stages 1 and 6 here), 
```r
dim_reduction_output$projection
```
|sample     |category |furthest_groups |   distance| normalised_distance|group |
|:----------|:--------|:---------------|----------:|-------------------:|:-----|
|SRR7610144 |appetite |s1-s6           | -2.8395655|           0.7192360|s2    |
|SRR7610145 |appetite |s1-s6           | -2.8303468|           0.7186562|s2    |
|SRR7610146 |appetite |s1-s6           | -2.8207618|           0.7180534|s3    |
|SRR7610147 |appetite |s1-s6           | -2.6803610|           0.7092227|s3    |
...
|SRR7610166 |vision   |s5-s1           | -4.127260|            0.1275470|s7    |
|SRR7610167 |vision   |s5-s1           | -3.394484|            0.1799185|s6    |
|SRR7610168 |vision   |s5-s1           | -5.911884|            0.0000000|s5    |
|SRR7610169 |vision   |s5-s1           | -3.772787|            0.1528812|s6    |
|SRR7610170 |vision   |s5-s1           | -4.073592|            0.1313827|s6    |
#
$information provides the category, number of variables (genes here) in the category also found in the counts table, the number of PCA dimensions retained and the variance they explain (default: threshold = 0.5 - so all sum_variance_kept_pcs > 0.5), and the correlation between the extracted value for each sample and the average gene expression of all genes in the category (default: spearman correlation test). 
```r
head(dim_reduction_output$information)
```
|category         | n_variables|method | num_pcs| sum_variance_kept_pcs| expr_pca_correlation| expr_pca_correlation_pvalue|
|:----------------|-----------:|:------|-------:|---------------------:|--------------------:|---------------------------:|
|appetite         |          80|PCA    |       3|             0.5578087|            0.6753247|                   0.0010605|
|digestion        |         121|PCA    |       2|             0.5502907|            0.9597403|                   0.0000051|
|gastrointestinal |          17|PCA    |       2|             0.6869954|            0.9883117|                   0.0000044|
|corticoids       |          28|PCA    |       3|             0.5676356|            0.3610390|                   0.1084280|
|thyroid          |          31|PCA    |       3|             0.5859959|            0.1116883|                   0.6287796|
|betaoxi          |          14|PCA    |       1|             0.6837123|            0.9324675|                   0.0000049|

*In this case, most biological categories show consistent expression across genes (_i.e.,_ genes show overall higher expression in a given sample). An exception is genes in the "thyroid" category here. Thyroid-related gene expression is often complex, and here different samples show different profiles with various genes acting in a complex heterogeneous pattern. Thus, for this biological category, the values extracted from **dim_reduction()** do not directly relate the overall gene expression level.*
#
$dimred_information provides key details for plotting PCA outputs for each category and visualising the main axis of variance (only for the first two PCs)
```r
head(dim_reduction_output$dimred_information)
```
|category         |method | num_pcs|       pc1|       pc2| centroid| maxvariancedirection|
|:----------------|:------|-------:|---------:|---------:|--------:|--------------------:|
|appetite         |PCA    |       3| 0.2365357| 0.1713829| -6.3e-17|            0.7167241|
|appetite         |PCA    |       3| 0.2365357| 0.1713829|  1.0e-15|           -0.6971527|
|appetite         |PCA    |       3| 0.2365357| 0.1713829| -1.2e-16|            0.0168694|
|digestion        |PCA    |       2| 0.4037782| 0.1465124|  3.2e-16|           -0.9881273|
|digestion        |PCA    |       2| 0.4037782| 0.1465124|  3.8e-16|           -0.1536375|
|gastrointestinal |PCA    |       2| 0.4675435| 0.2194519|  1.8e-16|           -0.9765832|

#
$pca provides the sample coordinates for the PCAs performed for each set of genes. Here is the result for the first 10 PCs of the appetite genes.
```r
head(dim_reduction_output$pca$appetite[,1"10])
```
|           |        PC1|      PC2|        PC3|        PC4|        PC5|        PC6|        PC7|        PC8|        PC9|       PC10|
|:----------|----------:|--------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
|SRR7610156 | -5.0416294| 3.781966|  1.4746764|  1.9508387|  0.1837337| -0.2140692| -0.4073693| -2.3561947|  1.2355917|  0.7849125|
|SRR7610157 | -7.0713019| 3.208837|  0.1010005|  0.3927382| -1.5593034| -2.0210512|  0.1321762|  1.0031456|  3.2866671| -1.2795689|
|SRR7610162 | -0.0442157| 3.485894|  1.1252763|  1.2029134|  1.2576794|  1.8854425| -0.7047708| -0.4350283| -1.5139857|  0.5091929|
|SRR7610144 | -2.0167031| 1.988757| -0.4551666| -0.5278917|  0.4099790|  1.2982684|  0.8531442|  0.0507996|  0.8293554| -2.5758929|
|SRR7610145 | -0.8675218| 3.211298|  1.7897778|  1.9636571|  1.7966160|  1.1435161|  1.6665057| -1.7197282|  1.2747393|  2.6163786|
|SRR7610163 | -2.3259845| 3.260200|  0.5103722|  1.7582506| -0.4610619|  0.4800184| -0.6355526| -0.0900312|  0.1237668| -1.4047689|

#
### Plotting
Using **plot_radar()** will plot the value extracted for each sample and biological category.
```r
plot_radar(data_input, dim_reduction_output)
```
![Radar plot for stage 1](example1/stage_1_raw.png)
*Figure 3: Radar plot for stage 1 of the _A. ocellaris_ dataset obtained with method = **"pca"**, without manually arranging the order of the biological categories around the radar.*

#
The order of the categories around the plot matches that of the list of biological categories provided to **import_data()**. This can be modified to change the order of the categories or remove some categories from the radar plot.
```r
category_list_names = c("glycolysis","betaoxi","krebs","appetite","cholesterol",
                        "fatty","digestion","gastrointestinal",
                        "corticoids","thyroid","ossification","vision")

category_list=cbind(category_list_names,c(1:length(category_list_names)))

colnames(category_list)=c("category","order")

category_list=as.data.frame(category_list)

plot_radar(data_input, dim_reduction_output, category_list = category_list)
```
![Radar plot for stage 1](example1/stage_1.png)
*Figure 4: Radar plot for stage 1 of the _A. ocellaris_ dataset obtained with method = **"pca"**.*

#
Users can also modify the size of the labels, control the order of each group, and display radar plots together.
```r
radars = plot_radar(data_input, dim_reduction_output, category_list, axis_label_size = 1, radar_label_size = 1)

ordered_list = c("s1","s2","s3","s4","s5","s6","s7")

wrap_plots(radars[ordered_list], ncol=4, nrow=2)
```
![Radar plots for each stage](example1/all_stages.png)
*Figure 5: Radar plot for all 7 stages of the _A. ocellaris_ dataset obtained with method = **"pca"**.*

### Dimensional reduction and plotting - scale


Below, we show the output using method = **"scale"**. For some biological categories, the results look similar to those on the radars obtained with method = **"pca"** - which is to be expected when the values extracted for each sample are closely linearly correlated to the mean expression in the sample (see $information output from __dim_reduction()__).
```r
dim_reduction_output = dim_reduction(
  data_input,
  method = "scale")

radars = plot_radar(data_input, dim_reduction_output, category_list, axis_label_size = 1, radar_label_size = 1)

wrap_plots(radars[ordered_list], ncol=4, nrow=2)
```
![Radar plots for each stage](example1/all_samples_scale.png)
*Figure 6: Radar plot for all 7 stages of the _A. ocellaris_ dataset obtained with method = **"scale"**.*

---
## Example #2 - method = "lda"
Here we go through an example pipeline using method = "lda".

We use the RNAseq data from the 3-time point developmental series under control (DMSO) and exposure treatments (2 substances at 2 concentrations each: 1.3 uM and 2.4 uM Sorafenib, 25 nM and 50 nM Rotenone) of zebrafish _Danio rerio_ (from [Nöth et al. (2025)](https://doi.org/10.1007/s00204-024-03944-7)). There are 45 samples in total (3 samples per stage x 3 stages x (1 control + 4 treatments)).

In this experiment, hour post fertilisation has a major footprint on gene expression profiles. Here, we run the pipeline with method = **"pca"** - and see that the values on the radar plot are strongly related to hour post fertilisation. Then, we use method = **"lda"** to extract the effect of the exposure treatments and their concentrations on various biological categories.

The biological categories are those defined for zebrafish in a manually curated database ([Herrera et al., 2025](https://doi.org/10.1002/jez.b.23299)].

```r
# load the package
library(RadarOmics)

# import data
data_input = import_data(expr_path = "vsd_zebrafish.csv", sample_meta_path = "sampleinfo_zebrafish.csv", gene_meta_path = "genelist_zebrafish.csv")
```
Expression data (or other tabular data), normalised for PCA use. For gene expression data, we recommend VSD normalisation with DESEq2.
```r
head(data_input$expr[,1:4])
```

|                   | 10_DMSO_Kontrolle2_96h_Eppi-12_Index-B5| 11_DMSO_Kontrolle3_96h_Eppi-13_Index-C5| 12_DMSO_Kontrolle4_96h_Eppi-14_Index-D5| 19_DMSO_Kontrolle3_36h_Eppi-28_Index-C6|
|:------------------|---------------------------------------:|---------------------------------------:|---------------------------------------:|---------------------------------------:|
|ENSDARG00000000001 |                                     145|                                     336|                                       0|                                     125|
|ENSDARG00000000002 |                                     719|                                    2285|                                       0|                                     456|
|ENSDARG00000000018 |                                     918|                                    2787|                                       0|                                    2997|
|ENSDARG00000000019 |                                    3901|                                   11648|                                       0|                                    4131|
|ENSDARG00000000068 |                                    1465|                                    4536|                                       0|                                     825|
|ENSDARG00000000069 |                                    2178|                                    6529|                                       0|                                    4625|

*Note that samples are columns, genes are rows.*

- Sample information with columns "sample", "group", and other information. One radar plot per "group" (unless other grouping is specified using _radar_grouping_ in **plot_radar()**) will be generated.
```r
head(data_input$sample_meta)
```
|sample                                   |group               | time_hpf|substance | concentration|substance_concentration |substance_hpe |
|:----------------------------------------|:-------------------|--------:|:---------|-------------:|:-----------------------|:-------------|
|51_DMSO_Kontrolle1_96h_Eppi-76_Index-D10 |96_DMSO_0_72        |       96|DMSO      |           0.0|DMSO_0                  |DMSO_96       |
|52_DMSO_Kontrolle2_96h_Eppi-77_Index-E10 |96_DMSO_0_72        |       96|DMSO      |           0.0|DMSO_0                  |DMSO_96       |
|53_DMSO_Kontrolle3_96h_Eppi-78_Index-F10 |96_DMSO_0_72        |       96|DMSO      |           0.0|DMSO_0                  |DMSO_96       |
|54_Sorafenib_EC50_96h_Eppi-83_Index-G10  |96_Sorafenib_1.3_72 |       96|Sorafenib |           1.3|Sorafenib_1.3           |Sorafenib_96  |
|55_Sorafenib_EC50_96h_Eppi-84_Index-H10  |96_Sorafenib_1.3_72 |       96|Sorafenib |           1.3|Sorafenib_1.3           |Sorafenib_96  |
|56_Sorafenib_EC50_96h_Eppi-85_Index-A11  |96_Sorafenib_1.3_72 |       96|Sorafenib |           1.3|Sorafenib_1.3           |Sorafenib_96  |

- Gene information with columns "gene" and "category". Users can manually reshuffle and filter categories before plotting the output of the package using the radar plot.
```r
data_input$gene_meta
```
|gene               |category          |
|:------------------|:-----------------|
|ENSDARG00000053449 |appetite          |
|ENSDARG00000035350 |appetite          |
|ENSDARG00000116735 |appetite          |
...
|ENSDARG00000075996 |brain_development |
|ENSDARG00000078701 |brain_development |
|ENSDARG00000018765 |brain_development |


### Dimensional reduction and plotting - PCA with "group"

Once the dataset is uploaded, we first explore the data with method = **"pca"** and plot the results.
```r
dim_reduction_output = dim_reduction(
  data_input,
  method = "pca")

radars=plot_radar(result, dim_reduction_output, category_list = category_list, axis_label_size=2.5, radar_label_size=3, radar_label_position = "top")

unique(result$sample_meta$group) # to control the order of the groups

wrap_plots(radars[unique(result$sample_meta$group)], ncol=5, nrow=3)
```
![Radar plot for all samples](example2/all_samples_pca.png)
*Figure 7: Radar plot for all groups (hour post fertilisation + substance + concentration) of the _D. rerio_ dataset obtained with method = **"pca"**, looking at position of samples along the main axis of variance running through all groups.*

On this plot, we see that hour post fertilisation (organised by row) exerts a major influence on biological processes, although some treatment-specific differences (organised by column) are already visible.

###
### Dimensional reduction and plotting - PCA with other grouping design

We can modify the function call by adding the argument _focus = "substance_concentration"_ which will use substance + concentration rather than groups to determine the main axis of variance when > 1 PC is retained.
```r
dim_reduction_output = dim_reduction(
  data_input,
  method = "pca",
  focus = "substance_concentration")

radars=plot_radar(result, dim_reduction_output, category_list = category_list, axis_label_size=2.5, radar_label_size=3, radar_label_position = "top")

unique(result$sample_meta$group) # to control the order of the groups

wrap_plots(radars[unique(result$sample_meta$group)], ncol=5, nrow=3)
```
![Radar plot for all samples](example2/all_samples_pca_substance_concentration.png)
*Figure 8: Radar plot for all groups (hour post fertilisation + substance + concentration) of the _D. rerio_ dataset obtained with method = **"pca"**, looking at position of samples along the main axis of variance running through all combinations of "substance_concentration".*

###
### Dimensional reduction and plotting - LDA

As an alternative to modifying the sample projection axes on the PCAs, in the case of complex and nested experimental treatments, the combination of PCA + LDA (method = "lda") can be used to better extract the footprint of a given treatment on expression profiles.

We thus test method = **"lda"** with _lda_focus = "substance_concentration"_, keeing the default _focus = "group"_. If the variance explained by PC1 > 80% of the total variance, the PC1 coordinate of each sample is retained. If not, a LDA is generated and the coordinates of each sample along the top LD dimensions representing at least 80% of the variance are used.

Using a high threshold for PCA variance allows to retain most of the inter-sample variance before feeding the coordinates to LDA. Users can test multiple thresholds and sets of _focus_ and _lda_focus_ and decide which are best suited for their dataset.

```r
dim_reduction_output = dim_reduction(
  data_input,
  method = "lda",
  lda_focus = "substance_concentration",
  threshold=0.8,
  lda_threshold = 0.8)

radars=plot_radar(result, dim_reduction_output, category_list = category_list, axis_label_size=2.5, radar_label_size=3, radar_label_position = "top")

unique(result$sample_meta$group) # to control the order of the groups

wrap_plots(radars[unique(result$sample_meta$group)], ncol=5, nrow=3)
```
![Radar plot for all samples](example2/all_samples_lda_08_08_substance_concentration.png)
*Figure 9: Radar plot for all groups (hour post fertilisation + substance + concentration) of the _D. rerio_ dataset obtained with method = **"lda"**, looking at position of samples along the main axis of variance running through all groups.*

In this instance, forcing the variance to capture the footprint of the substance + concentration allows to better discern to what extent biological processes are modified by each treatment, across the different time points.

Note that the values extracted by **dim_reduction()** do not necessarily reflect the overall gene expression level (_i.e.,_ higher _vs._ lower) but rather show how distant samples are from each other in the multidimensional space.

We recommend checking the various outputs of **dim_reduction()** to fully inspect the results and interpret the radar plots.

---

### Alternative plot grouping

Note that users can customise the grouping of the samples on radar plots using the argument _radar_grouping_ in **plot_radar()**
```r
radars=plot_radar(result, dim_reduction_output,radar_grouping = "substance_concentration",
category_list = category_list,radar_label_size=3,axis_label_size=2.5,radar_label_position = "top")

wrap_plots(radars[unique(result$sample_meta$substance_concentration)],ncol=5,nrow=1)
```
![radar based on substance + concentration grouping rather than all groups](example2/substance_concentration_radar.png)
*Figure 10: Radar plot for all combinations of + substance + concentration of the _D. rerio_ dataset obtained with method = **"lda"**, looking at position of samples along the main axis of variance running through all groups.*

---
## Data inspection

### Basic information
Basic information about the number of genes/proteins/others used in each category and statistics about the values extracted by **dim_reduction()** are provided by the outputs $projection, $information, and $dimred_information. See examples of output in the section [Example #1 | Dimensional reduction - PCA](#dimensional-reduction---PCA).

### PCA and LDA plots
In addition to yielding simple markers such as the correlation between the PCA/LDA-derived values for each sample and the expression level in each category in $information from the output of **dim_reduction()**, we also provide options to plot the PCAs and LDAs and visualise the main axis of variance when > 1 dimension is retained (if the first dimension passes the user-defined variance threshold, the PC1 or LD1 coordinate is used to derive values for each sample and biological category).

Still using the dataset from zebrafish, we use **plot_dimension()** to visualise the first two dimensions of the dimensional reduction approach chosen (here, PCA if the variance explained by PC1 is >0.8, if not, LDA) and the main axis of variance (when more than 1 LD is retained, if not the LD1 coordinate is used by **dim_reduction()** to extract the values to plot on the radars for each sample and biological category).

```r
plots = plot_dimensions(data_input, dim_reduction_output,
                      colour="substance_hpf", shape="hpf",
                      point_size=3, colour_palette = c("skyblue1","skyblue2","skyblue3","orange","orange3","orange4","pink1","pink2","red3"))

plots$vision # PC1 met the threshold
plots$glycolysis # LD1 met the threshold
plots$cholesterol # LD1 met the threshold
plots$thyroid # LD1 + 2 met the threshold, thus necessitating the calculation of the main axis of variance across groups of samples.
```
Find examples of plots and some comments below.


![PC1 and 2 for vision](example2/pca_vision.png)
*Figure 11: PC1 and 2 obtained from the vision-related genes.*

Here, we can see the signature of hour post fertilisation on visual genes, and discern an effect of Sorafenib on vision genes at 96 hours.

###
![LD1 and 2 for glycolysis](example2/lda_glycolysis.png)
*Figure 12: LD1 and 2 obtained from the glycolysis-related genes, with _lda_focus = "substance_concentration"_.*

Here, we can see the value of moving from PCA to LDA to minimise the effect of hour post fertilisation and better isolate that of treatments on glycolysis-related genes.

###
![LD1 and 2 for cholesterol](example2/lda_cholesterol.png)
*Figure 13: LD1 and 2 obtained from the cholesterol synthesis-related genes, with _lda_focus = "substance_concentration"_.*

Here, we can see the value of moving from PCA to LDA to minimise the effect of hour post fertilisation on profiling and better isolate that of treatments, particularly the effect of Sorafenib on cholesterol synthesis-related genes at 36 and 96 hours.

###
![LD1 and 2 for thyroid](example2/lda_thyroid.png)
*Figure 14: LD1 and 2 obtained from the thyroid-related genes, with _lda_focus = "substance_concentration"_.*

This is a case where using _focus = "group"_ leads the main axis of variance to represent hour post fertilisation more than treatment, as the groups are defined based on the combination of hour post fertilisation + treatment + concentration. Using a different _focus_ that does not include hour post fertilisation may lead to better disentangling of the effect of treatment.

## Boxplots of values yielded by **dim_reduction()** and statistical testing
The values extracted for each sample can be plotted on boxplots with  **plot_boxplot()**, which also displays statistical test results (aa, ab, bc, etc.). These boxplots can be used as criteria for refining the choice of _focus_, _lda_focus_, _threshold_, and _lda_threshold_ to use as part of **dim_reduction()**.

