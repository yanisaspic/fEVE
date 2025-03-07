"Functions called to set-up the data of a scEVE clustering recursion.

	2025/03/06 @yanisaspic"

sceve_GetSelectedFeatures <- function(dataset, params, n_features=1000) {
  #' Get the n most variable genes in a scRNA-seq dataset.
  #'
  #' The variable genes are identified by calling the function FindVariableFeatures() of Seurat.
  #'
  #' @param dataset a scRNA-seq dataset of raw count dataset, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param n_features the number of highly variable genes to sample.
  #'
  #' @return a vector of genes.
  #'
  #' @import Seurat
  #'
  #' @export
  #'
  SeurObj <- Seurat::CreateSeuratObject(dataset)
  SeurObj <- Seurat::FindVariableFeatures(SeurObj, nfeatures = n_features)
  return(Seurat::VariableFeatures(SeurObj))
}