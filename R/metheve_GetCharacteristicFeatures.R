"Functions called to predict marker methylation sites.

	2025/04/22 @yanisaspic"

get_differential_methylation <- function(site, filtered_dataset, groups) {
  #' Get the log2 fold change, as well as the FDR, of a methylation site.
  #'
  #' @param site a character.
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #' @param groups an ordered vector of integers.
  #'
  #' @return a data.frame with three columns: `site`, `logFC` and `pval`.
  #'
  #' @import stats
  #'
  in_methylation <- as.numeric(filtered_dataset[site, groups])
  out_methylation <- as.numeric(filtered_dataset[site, !groups])
  logfc <- log2(mean(in_methylation) / mean(out_methylation))
  pval <- stats::wilcox.test(in_methylation, out_methylation)$p.value
  differential_methylation <- data.frame(site=site, logFC=logfc, pval=pval)
  return(differential_methylation)
}

metheve_GetCharacteristicFeatures <- function(cluster, selected_data, params,
                                            logFC_threshold=2, pvalue_threshold=0.001) {
  #' Get differentially methylated features by using the Wilcoxon test.
  #' Only sites with abs(logFC) > 2 and FDR < 0.001 are returned.
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
  is_in_cluster <- function(sample) {ifelse(sample %in% cluster$samples, TRUE, FALSE)}
  groups <- is_in_cluster(samples_of_recursion)

  f <- function(site) {get_differential_methylation(site, selected_data$dataset, groups)}
  tmp <- lapply(X=rownames(selected_data$dataset), FUN=f)
  tmp <- do.call(rbind, tmp)
  tmp$FDR <- stats::p.adjust(tmp$pval, method="BH")
  tmp <- tmp[(abs(tmp$logFC)>logFC_threshold) & (tmp$FDR<pvalue_threshold),]
  markers <- stats::setNames(tmp$logFC, tmp$site)
  markers <- markers[order(-markers)]
  return(markers)
}
