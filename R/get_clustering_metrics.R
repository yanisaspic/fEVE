"Miscellaneous functions called by multiple scripts in the fEVE package.

	2025/05/09 @yanisaspic"

get_data_bluster <- function(dataset) {
  #' Get a matrix usable for the functions of bluster.
  #'
  #' @param dataset a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #'
  #' @return an object associating samples and selected features to their reduced dimensions.
  #'
  #' @import scater
  #' @import scran
  #' @import scuttle
  #' @import SingleCellExperiment
  #'
  data <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(dataset)))
  data <- scuttle::logNormCounts(data)
  variable_genes <- scran::getTopHVGs(scran::modelGeneVar(data), n=5000)
  set.seed(1)
  data <- scater::runPCA(data, ncomponents=20, subset_row=variable_genes)
  data <- SingleCellExperiment::reducedDim(data)
  return(data)
}

get_intrinsic_clustering_metrics <- function(dataset, preds) {
  #' Using intrinsic metrics, measure a clustering performance.
  #'
  #' @param dataset a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param preds a named vector associating cells to their predicted clusters.
  #'
  #' @return a named vector of with two names: `nPurity` and `SI`.
  #'
  #' @import bluster
  #'
  n_clusters_predicted <- length(unique(preds))
  if (n_clusters_predicted < 2) {return(c("nPurity"=NA, "SI"=NA))}
  data_bluster <- get_data_bluster(dataset)
  neighborhood_purity <- bluster::neighborPurity(data_bluster, preds)
  silhouette_index <- bluster::approxSilhouette(data_bluster, preds)
  intrinsic_clustering_metrics <- c("nPurity"=mean(neighborhood_purity$purity),
                                    "SI"=mean(silhouette_index$width))
  return(intrinsic_clustering_metrics)
}

get_clustering_metrics <- function(data, preds) {
  #' Using both intrinsic and extrinsic clustering metrics, measure a clustering performance.
  #'
  #' Extrinsic metrics compare cluster predictions to the cell annotations of the dataset.
  #' Intrinsic metrics compare the gene expression of cells in and out of their clusters.
  #' The `ARI` and the `NMI` are extrinsic metrics.
  #' The `nPurity` and the `SI` are intrinsic metrics.
  #' For every metric, higher is better and the maximum value is 1.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a dataset, without selected features. Its rows are features and its columns are samples.
  #' `ground_truth` is a named factor associating samples to their cluster annotations.
  #' @param preds a named vector associating cells to their predicted clusters.
  #'
  #' @return a named vector with four names: `ARI`, `NMI`, `nPurity` and `SI`.
  #'
  #' @import aricode
  #'
  if (length(preds) < 2) {return(c("ARI"=NA, "NMI"=NA, "nPurity"=NA, "SI"=NA))}
  ground_truth <- data$ground_truth[names(preds)]
  dataset <- data$dataset[, names(preds)]
  clustering_metrics <- c("ARI"=aricode::ARI(ground_truth, preds),
                          "NMI"=aricode::NMI(ground_truth, preds),
                          get_intrinsic_clustering_metrics(dataset, preds))
  return(clustering_metrics)
}

get_clustering_metrics_trifecta <- function(data, preds, method_label) {
  #' Using both intrinsic and extrinsic clustering metrics, measure the clustering performances
  #' of a fEVE clustering (with and without its leftover clusters), and the ground truth.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a dataset, without selected features. Its rows are features and its columns are samples.
  #' `ground_truth` is a named factor associating samples to their cluster annotations.
  #' @param preds a named vector associating cells to their predicted clusters.
  #' @param method_label a character.
  #'
  #' @return a data.frame with six columns: `method`, `ARI`, `NMI`, `nPurity`, `SI` and `n_samples`.
  #'
  #' @export
  #'
  population_is_leftover <- function(population) {
    substr(population, nchar(population), nchar(population))=="L"}
  leftover_populations <- Filter(population_is_leftover, unique(preds))
  robust_populations <- preds[!preds %in% leftover_populations]

  x <- c("method"=method_label, get_clustering_metrics(data, preds), "n_samples"=length(preds))
  y <- c("method"=glue::glue("{method_label}*"), get_clustering_metrics(data, robust_populations), "n_samples"=length(robust_populations))
  z <- c("method"="ground_truth", get_clustering_metrics(data, data$ground_truth), "n_samples"=length(data$ground_truth))
  benchmark <- as.data.frame(rbind(x,y,z))
  return(benchmark)
}