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
Here is an example of usage based on RNAseq data from the 7-stage developmental series of the false clownfish _Amphiprion ocellaris_ (from [Roux et al. (2023)](https://doi.org/10.1016/j.celrep.2023.112661).

We use radar plots to summarise the gene expression profile of each sample at each stage for a pre-defined set of biological processes.

We are supplying:
- a variance-stabilisation transformed count table (vsd) obtained using DESEq2 with genes as rows, samples as columns
- sample information with two columns: samples and their grouping (here, developmental stage, from stage 1 to stage 7)
- gene list with two columns: genes and their categories (here, various biological processes, _e.g.,_ glycolysis, Krebs cycle, phototransduction)

1) load the package
```r
library(RadarOmics)
```
2) Upload the correctly-formatted required input files
```r

```
