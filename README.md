# RadarOmics

RadarOmics is an R package for summarising **'omics datasets** through **dimensional reduction**  and visualising the output across samples, treatments, and biological processes with **radar plots**.
It is designed to handle **complex experimental designs** with >2 treatments, and can handle nested designs (_e.g._, multiple developmental stages and chemical exposure treatments). It works with gene expression, protein expression, or similar tabular datasets.
#
---

# Installation

You can install the development version directly from GitHub:

```r
# Install remotes if not already installed
install.packages("remotes")

# Install RadarOmics from GitHub
remotes::install_github("Emma-Gairin/RadarOmics", auth_token = "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")
```

### Implementation and example of use

Here we demonstrate the use of **RadarOmics** to summarise the gene expression profile of each sample within different groups for a pre-defined set of biological processes.
We use as an example the RNAseq data from the 7-stage developmental series of the false clownfish _Amphiprion ocellaris_ (from [Roux et al. (2023)](https://doi.org/10.1016/j.celrep.2023.112661)).

#### Load the package
```r
# load the package
library(RadarOmics)
```

#### Import the data.
We are supplying:
- a variance-stabilisation transformed count table (vsd) obtained using DESEq2 with genes as rows, samples as columns
- sample information with two columns: samples and their grouping (here, developmental stage, from stage 1 to stage 7)
- gene list with two columns: genes and their categories (here, various biological processes, _e.g.,_ glycolysis, Krebs cycle, phototransduction)

```r
# import data
data_input = import_data(expr_path = "vsd_ocellaris.csv", sample_meta_path = "sampleinfo_ocellaris.csv", gene_meta_path = "genelist_ocellaris.csv")
```
Here is what the data should look like.
- Expression data (or other tabular data), normalised for PCA use. For gene expression data, we recommend VSD normalisation with DESEq2.
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

- Sample information with columns "sample" and "group". One radar plot per "group" will be generated.
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

- Gene information with columns "gene" and "category". Users can manually reshuffle and filter categories before plotting the output of the package using the radar plot.
```r
head(data_input$gene_meta)
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

#### Dimensional reduction
Once the dataset is uploaded, we can run the PCA and extract reduced coordinates from each sample and each biological category based on top PC dimensions representing e.g., 40 % of variance (defined by threshold = 0.4).
```r
dim_reduction_output = dim_reduction(
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

