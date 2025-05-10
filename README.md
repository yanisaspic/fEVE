# fEVE: An Alternative Framework to -Omics Ensemble Clustering

fEVE is an alternative approach to integrate multiple clustering results of -omics datasets. Instead of minimizing the differences across clustering predictions, fEVE describes these differences, and it leverages them to identify clusters that are robust to the method used.

Currently, two instances of the fEVE framework have been implemented:

-   **scEVE:** a specialized instance that analyzes single-cell transcriptomics datasets, using four scRNA-seq clustering methods. An in-depth justification of our methodological choices is available in the scEVE Method paper (in press).

-   **fEVE:** a generic instance that analyzes -omics datasets, using four data-agnostic clustering methods. A brief overview of our methodological choices is available in the fEVE Application Notes paper (in writing).

## Installation

You can install fEVE from Github with:

```{r}
install.packages("devtools")
devtools::install_github("yanisaspic/fEVE", dependencies=TRUE)
```

## Overview of the fEVE framework

-   A complete overview of the recursive clustering framework is available in the vignette `feve.Rmd`.

-   An overview of the on-demand scRNA-seq benchmark is available in the vignette `get_benchmark.Rmd`.

-   An overview of the analyses' reports generated with fEVE is available in the vignette `get_analysis_report.Rmd`
