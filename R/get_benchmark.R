"Functions called to benchmark an instance of the fEVE framework.

	2025/03/06 @yanisaspic"

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
  #' @param preds a named factor associating cells to their predicted clusters.
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
  #' @param preds a named factor associating cells to their predicted clusters.
  #'
  #' @return a named vector with four names: `ARI`, `NMI`, `nPurity` and `SI`.
  #'
  #' @import aricode
  #'
  #' @export
  #'
  if (length(preds) < 2) {return(c("ARI"=NA, "NMI"=NA, "nPurity"=NA, "SI"=NA))}
  ground_truth <- data$ground_truth[names(preds)]
  dataset <- data$dataset[, names(preds)]
  clustering_metrics <- c("ARI"=aricode::ARI(ground_truth, preds),
                          "NMI"=aricode::NMI(ground_truth, preds),
                          get_intrinsic_clustering_metrics(dataset, preds))
  return(clustering_metrics)
}

get_benchmark <- function(data, params, method_label) {
  #' Using computational as well as intrinsic and extrinsic clustering metrics, measure
  #' the performance of an instance of the fEVE framework on a dataset.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a dataset, without selected features. Its rows are features and its columns are samples.
  #' `ground_truth` is a named factor associating samples to their cluster annotations.
  #' @param params a list of parameters (cf. `feve::get_default_parameters()`).
  #' @param method_label a character.
  #'
  #' @return a data.frame with eight columns: `method`, `time (s)`,
  #' `peak_memory_usage (Mb)`, `ARI`, `NMI`, `nPurity`, `SI` and `n_samples`.
  #'
  #' @import glue
  #'
  #' @export
  #'
  get_memory_usage <- function(memory) {memory[[11]] + memory[[12]]}
  memory_usage_init <- get_memory_usage(gc(reset=TRUE))
  time_init <- Sys.time()
  results <- feve(data$dataset, params, figures=FALSE, sheets=FALSE)
  time <- as.numeric(Sys.time() - time_init, units="secs")
  peak_memory_usage <- get_memory_usage(gc()) - memory_usage_init

  population_is_leftover <- function(population) {
    (results$records$meta[population, "robustness"] == 0)}
  leftover_populations <- Filter(population_is_leftover, levels(results$preds))

  benchmark_all <- c("method"=method_label, "time (s)"=time,
                     "peak_memory_usage (Mb)"=peak_memory_usage,
                     get_clustering_metrics(data, results$preds),
                     "n_samples"=length(results$preds))
  benchmark_robust <- c("method"=glue::glue("{method_label}*"), "time (s)"=time,
                        "peak_memory_usage (Mb)"=peak_memory_usage,
                        get_clustering_metrics(data, results$preds[!results$preds %in% leftover_populations]),
                        "n_samples"=length(results$preds[!results$preds %in% leftover_populations]))
  ground_truth <- c("method"="ground_truth", "time (s)"=NA, "peak_memory_usage (Mb)"=NA,
                    get_clustering_metrics(data, data$ground_truth),
                    "n_samples"=length(data$ground_truth))
  benchmark <- as.data.frame(rbind(benchmark_all, benchmark_robust, ground_truth))
  return(benchmark)
}
