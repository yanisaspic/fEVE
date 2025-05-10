"Functions called to predict differentially expressed features.

	2025/05/09 @yanisaspic"

get_differential_feature <- function(feature, filtered_dataset, groups) {
  #' Get the log2 fold change, as well as the FDR, of a feature.
  #'
  #' @param feature a character.
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #' @param groups an ordered vector of integers.
  #'
  #' @return a data.frame with three columns: `feature`, `logFC` and `pval`.
  #'
  #' @import stats
  #'
  ingroup_feature <- as.numeric(filtered_dataset[feature, groups])
  outgroup_feature <- as.numeric(filtered_dataset[feature, !groups])
  logfc <- log2(mean(ingroup_feature) / mean(outgroup_feature))
  pval <- stats::wilcox.test(ingroup_feature, outgroup_feature)$p.value
  result <- data.frame(feaure=feature, logFC=logfc, pval=pval)
  return(result)
}

feve_GetCharacteristicFeatures <- function(cluster, selected_data, params,
                                           logFC_threshold=2, pvalue_threshold=0.001) {
  #' Get differentially expressed features by using the Wilcoxon test.
  #' Only sites with abs(logFC) > 2 and FDR < 0.001 are returned.
  #'
  #' @param cluster a named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param logFC_threshold a numeric.
  #' @param pvalue_threshold a numeric.
  #'
  #' @return a named vector associating differentially expressed features to their log2 fold-changes.
  #'
  #' @import stats
  #'
  #' @export
  #'
  samples_of_recursion <- colnames(selected_data$dataset)
  is_in_cluster <- function(sample) {ifelse(sample %in% cluster$samples, TRUE, FALSE)}
  groups <- is_in_cluster(samples_of_recursion)

  f <- function(feature) {get_differential_feature(feature, selected_data$dataset, groups)}
  tmp <- lapply(X=rownames(selected_data$dataset), FUN=f)
  tmp <- do.call(rbind, tmp)
  tmp$FDR <- stats::p.adjust(tmp$pval, method="BH")
  tmp <- tmp[(abs(tmp$logFC)>logFC_threshold) & (tmp$FDR<pvalue_threshold),]
  markers <- stats::setNames(tmp$logFC, tmp$feature)
  markers <- markers[order(-markers)]
  return(markers)
}
