# fEVE: An Alternative Framework to -Omics Ensemble Clustering

fEVE is an alternative approach to integrate multiple clustering results, generated with -omics datasets. 
Instead of minimizing the differences of multiple clustering results generated with a panel of methods fEVE describes these differences, and it leverages them to identify clusters that are robust to the method used.

Currently, one instance of the fEVE framework has been implemented:
    - scEVE: an instance that analyzes single-cell transcriptomics datasets, using four clustering methods.

fEVE is maintained by Asloudj Yanis [yanis.asloudj@u-bordeaux.fr].

## Installation

You can install fEVE from Github with:

```{r}
install.packages("devtools")
devtools::install_github("yanisaspic/fEVE", dependencies=TRUE)
```

## Overview of the fEVE framework

A complete overview of the framework is available in its associated vignette.