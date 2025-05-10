"Function called to benchmark an instance of the fEVE framework.

	2025/05/09 @yanisaspic"

get_benchmark <- function(data, params, method_label="fEVE") {
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
  computational_metrics <- c("time (s)" = as.numeric(Sys.time() - time_init, units="secs"),
                             "peak_memory_usage (Mb)" = get_memory_usage(gc()) - memory_usage_init)

  benchmark <- get_clustering_metrics_trifecta(data, results$preds, method_label)
  for (cm in names(computational_metrics)) {
    benchmark[, cm] <- computational_metrics[cm]
    benchmark[benchmark$method=="ground_truth", cm] <- NA}
  return(benchmark)
}