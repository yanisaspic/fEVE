"Functions called to predict marker features.

	2025/03/05 @yanisaspic"

get_differential_expression <- function(filtered_dataset, groups) {
  #' Get the features differentially expressed in an in-group (1),
  #' with regards to an out-group (0).
  #'
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #' @param groups an ordered vector of integers.
  #'
  #' @return a data.frame associating genes to their log2 fold change and FDR.
  #'
  #' @import edgeR
  #'
  tmp <- edgeR::DGEList(counts=filtered_dataset, group=groups)
  tmp <- edgeR::estimateDisp(tmp)
  test_results <- edgeR::exactTest(tmp)
  test_results <- edgeR::topTags(test_results, n=nrow(test_results$table))
  test_results <- test_results$table
  test_results <- test_results[order(test_results$logFC, decreasing=TRUE), ]
  return(test_results)
}

breve_GetCharacteristicFeatures <- function(cluster, selected_data, params,
                                            logFC_threshold=2, pvalue_threshold=0.001) {
  #' Get differentially expressed features by calling the function `FindMarkers` of the edgeR package.
  #' Only genes with abs(logFC) > 2 and FDR < 0.001 are returned.
  #'
  #' @param cluster a named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param logFC_threshold a numeric.
  #' @param pvalue_threshold a numeric.
  #'
  #' @return a named vector associating over-expressed features to their log2 fold-changes.
  #'
  #' @import stats
  #'
  #' @export
  #'
  samples_of_recursion <- colnames(selected_data$dataset)
  is_in_cluster <- function(sample) {ifelse(sample %in% cluster$samples, 1, 0)}
  n_samples_in <- sum(is_in_cluster(samples_of_recursion))
  n_samples_out <- length(samples_of_recursion) - n_samples_in
  if (n_samples_in < 3 | n_samples_out < 3) {return(c())}
  groups <- factor(is_in_cluster(samples_of_recursion))

  tmp <- get_differential_expression(selected_data$dataset, groups)
  tmp <- tmp[(abs(tmp$logFC)>logFC_threshold) & (tmp$FDR<pvalue_threshold),]
  markers <- stats::setNames(tmp$logFC, rownames(tmp))
  markers <- markers[order(-markers)]
  return(markers)
}
