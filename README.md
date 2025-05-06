# fEVE: An Alternative Framework to -Omics Ensemble Clustering

fEVE is an alternative approach to integrate multiple clustering results, generated with -omics datasets. Instead of minimizing the differences of multiple clustering results generated with a panel of methods fEVE describes these differences, and it leverages them to identify clusters that are robust to the method used.

Currently, two instances of the fEVE framework have been implemented:

-   **scEVE:** an instance that analyzes single-cell transcriptomics datasets, using four clustering methods. An in-depth justification of our methodological choices is available in the scEVE Method paper (in press).

-   **brEVE:** an instance that analyzes bulk transcriptomics datasets, using four clustering methods. A brief overview of our methodological choices is available in the fEVE Application Notes paper (in writing).

## Installation

You can install fEVE from Github with:

```{r}
install.packages("devtools")
devtools::install_github("yanisaspic/fEVE", dependencies=TRUE)
```

## Overview of the fEVE framework

A complete overview of the framework is available in the vignette `feve.Rmd`.\
An overview of the on-demand scRNA-seq benchmark is also available in the vignette `benchmark.Rmd`.
