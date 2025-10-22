# RadarOmics

RadarOmics is an R package for summarising **'omics datasets** through **dimensional reduction**  and visualising the output across samples, treatments, and biological processes with **radar plots**.
It is designed to handle **complex experimental designs** with >2 treatments, and can handle nested designs (_e.g._, multiple developmental stages and chemical exposure treatments). It works with gene expression, protein expression, or similar tabular datasets.
#
---

## Installation

You can install the development version directly from GitHub:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install RadarOmics from GitHub
remotes::install_github("Emma-Gairin/RadarOmics", auth_token = "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")
```

## Implementation
1) Load the package
```r
library(RadarOmics)
```
2) Upload the correctly-formatted required input files
```r

```


## Example of use
Here is an example based on RNAseq data from the 7-stage developmental series of the false clownfish _Amphiprion ocellaris_ (from [Roux et al. (2023)](https://doi.org/10.1016/j.celrep.2023.112661)).

We use radar plots to summarise the gene expression profile of each sample at each stage for a pre-defined set of biological processes.

We are supplying:
- a variance-stabilisation transformed count table (vsd) obtained using DESEq2 with genes as rows, samples as columns
- sample information with two columns: samples and their grouping (here, developmental stage, from stage 1 to stage 7)
- gene list with two columns: genes and their categories (here, various biological processes, _e.g.,_ glycolysis, Krebs cycle, phototransduction)

```r
# load the package
library(RadarOmics)

# import data
data_input <- import_data(expr_path = "vsd_ocellaris.csv", sample_meta_path = "sampleinfo_ocellaris.csv", gene_meta_path = "genelist_ocellaris.csv")

# run PCA and extract reduced coordinates from each sample and each biological category based on top PC dimensions representing e.g., 40 % of variance (defined by threshold = 0.4)

dim_reduction_output=dim_reduction(
  data_input,
  method = "pca",
  threshold = 0.4, focus= "group") # group here is the developmental stage
```

**dim_reduction()** yields multiple objects.
By using scaling, PCA, or LDA, it obtains a value between 0 and 1 for each sample.
$projection provides, for each sample and biological category, the two furthest groups (e.g., stages 1 and 6 here), 
```r
head(dim_reduction_output$projection)
```
|sample     |category |furthest_groups |   distance| normalised_distance|group |
|:----------|:--------|:---------------|----------:|-------------------:|:-----|
|SRR7610144 |appetite |s6-s1           |  2.8323352|           0.7202689|s2    |
|SRR7610145 |appetite |s6-s1           |  2.8659379|           0.7223736|s2    |
|SRR7610146 |appetite |s6-s1           |  2.8315239|           0.7202181|s3    |
|SRR7610147 |appetite |s6-s1           |  2.7019725|           0.7121038|s3    |
|SRR7610148 |appetite |s6-s1           |  0.7051692|           0.5870365|s3    |
|SRR7610149 |appetite |s6-s1           | -0.2501827|           0.5271992|s4    |

$information provides the number of dimensions retained from the PCA (and LDA if using LDA), the variance described by PC1 and PC2, the number of genes in the category, and the correlation between the extracted value for each sample and the average gene expression of all genes in the category (default: spearman). "flipped" is TRUE is the normalised_distance reported in $projection was reversed (i.e., if 
```r
head(dim_reduction_output$information)
```

|category         |method | num_pcs| sum_variance_kept_pcs|       pc1|       pc2| n_genes| expr_pca_correlation| expr_pca_correlation_pvalue|
|:----------------|:------|-------:|---------------------:|---------:|---------:|-------:|--------------------:|---------------------------:|
|appetite         |PCA    |       2|             0.4079186| 0.2365357| 0.1713829|      80|            0.4194805|                   0.0596069|
|digestion        |PCA    |       1|             0.4037782| 0.4037782| 0.1465124|     121|            0.9883117|                   0.0000044|
|gastrointestinal |PCA    |       1|             0.4675435| 0.4675435| 0.2194519|      17|            0.9896104|                   0.0000044|
|corticoids       |PCA    |       2|             0.4135218| 0.2340784| 0.1794435|      28|            0.3623377|                   0.1070934|
|thyroid          |PCA    |       2|             0.4327358| 0.2307809| 0.2019549|      31|            0.5961039|                   0.0051323|
|betaoxi          |PCA    |       1|             0.6837123| 0.6837123| 0.1197882|      14|            0.9337662|                   0.0000050|

$pca_information provides key details for plotting PCA outputs for each category and visualising the main axis of variance (only for the first two PCs)
```r
head(dim_reduction_output$pca_information
```
|category         |method | num_pcs|       pc1|       pc2| centroid| maxvariancedirection|
|:----------------|:------|-------:|---------:|---------:|--------:|--------------------:|
|appetite         |PCA    |       2| 0.2365357| 0.1713829|        0|           -0.7147671|
|appetite         |PCA    |       2| 0.2365357| 0.1713829|        0|            0.6993625|
|digestion        |PCA    |       1| 0.4037782| 0.1465124|       NA|                   NA|
|gastrointestinal |PCA    |       1| 0.4675435| 0.2194519|       NA|                   NA|
|corticoids       |PCA    |       2| 0.2340784| 0.1794435|        0|            0.0238073|
|corticoids       |PCA    |       2| 0.2340784| 0.1794435|        0|            0.9997166|
