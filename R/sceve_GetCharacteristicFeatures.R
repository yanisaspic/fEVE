"Functions called to predict marker genes.

	2025/03/05 @yanisaspic"

sceve_GetCharacteristicFeatures <- function(cluster, selected_data, params, FC_threshold=4, pvalue_threshold=0.001) {
  #' Get marker genes by calling the function `FindMarkers` of the Seurat package.
  #' Only genes with log2FC > 4 and p-values (corrected with Bonferonni) < 0.001 are returned.
  #'
  #' @param cluster a named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param FC_threshold a numeric.
  #' @param pvalue_threshold a numeric.
  #'
  #' @return a named vector associating marker genes to their log2 fold-changes.
  #'
  #' @import Seurat
  #' @import stats
  #'
  samples_of_recursion <- colnames(selected_data$dataset)
  is_in_cluster <- function(cell) {ifelse(cell %in% cluster$samples, 1, 0)}
  n_samples_in <- sum(is_in_cluster(samples_of_recursion))
  n_samples_out <- length(samples_of_recursion) - n_samples_in
  if (n_samples_in < 3 | n_samples_out < 3) {return(c())}

  Seurat::Idents(object=selected_data$SeuratObject) <- factor(is_in_cluster(samples_of_recursion))
  markers <- Seurat::FindMarkers(selected_data$SeuratObject, ident.1=1)
  markers <- markers[(markers$avg_log2FC > FC_threshold) & (markers$p_val_adj < pvalue_threshold), ]
  markers <- stats::setNames(markers$avg_log2FC, rownames(markers))
  markers <- markers[order(-markers)]
  return(markers)
}