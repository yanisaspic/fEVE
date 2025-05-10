"Main functions called to condut a fEVE clustering analysis.

	2025/05/09 @yanisaspic"

feve_recursion <- function(population, dataset_init, SeuratObject_init, records, params, figures, sheets) {
  #' Attempt to cluster a specific population using the fEVE algorithm.
  #'
  #' @param population a character. It corresponds to the population that scEVE will attempt to cluster.
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param SeuratObject_init a SeuratObject generated from dataset_init, on which
  #' the function RunUMAP() of Seurat has been applied already.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #' @param sheets a boolean that indicates if the results of the clustering recursion should be saved in Excel sheets.
  #'
  #' @return a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @import openxlsx
  #'
  while(TRUE) {
    selected_data <- get_selected_data(population, dataset_init, SeuratObject_init, records, params, figures)
    if (length(selected_data) == 0) {break()} # the population is too small, and an empty list is returned.
    base_clusters <- get_base_clusters(population, selected_data, params, figures)
    robust_clusters <- get_robust_clusters(population, base_clusters, selected_data, records, params, figures)
    if (length(robust_clusters) == 0) {break()} # robust clusters are not predicted, and an empty list is returned.
    characterized_clusters <- get_characterized_clusters(population, robust_clusters, selected_data, params, figures)
    if (length(characterized_clusters) == 0) {break()}  # the clusters are homogenous, and an empty list is returned.
    records <- get_recorded_recursion(population, characterized_clusters, selected_data, records, params)
    break()}
  records$meta[population, "clustering_status"] <- "COMPLETE"
  if (sheets) {openxlsx::write.xlsx(records, params$sheets_path, rowNames=TRUE)}
  if (figures) {merge_drawings(population, params)}
  return(records)
}

feve_main <- function(population, dataset_init, records, params, figures, sheets) {
  #' Conduct a clustering analysis with the fEVE framework, starting from a target population.
  #'
  #' @param population a character. It corresponds to the population that scEVE will attempt to cluster.
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #' @param sheets a boolean that indicates if the results of the clustering recursion should be saved in Excel sheets.
  #'
  #' @return a named list, with two elements: `records` and `preds`.
  #' `records` is a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' `preds` is a named vector associating samples to their predicted clusters.
  #'
  #' @import openxlsx
  #'
  if (figures) {
    dir.create(params$figures_path)
    SeuratObject_init <- get_SeuratObject_init(dataset_init)}
  else {SeuratObject_init <- NA}

  while (!is.na(population)) {
    records <- feve_recursion(population, dataset_init, SeuratObject_init, records, params, figures, sheets)
    population <- get_pending_population(records$meta)}

  feature_is_characteristic <- function(feature) {sum(abs(feature)) > 0}
  records$features <- records$features[apply(X=records$features, MARGIN=1, FUN=feature_is_characteristic),]
  if (sheets) {openxlsx::write.xlsx(records, params$sheets_path, rowNames=TRUE)}

  results <- list(records=records, preds=get_leaf_clusters(records$samples))
  return(results)
}

feve <- function(dataset_init, params, figures=TRUE, sheets=TRUE, init_population="C") {
  #' Conduct a clustering analysis with the fEVE framework, starting from the initial population C.
  #'
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #' @param sheets a boolean that indicates if the results of the clustering recursion should be saved in Excel sheets.
  #' @param init_population a character (without `.`).
  #'
  #' @return a named list, with two elements: `records` and `preds`.
  #' `records` is a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' `preds` is a named vector associating samples to their predicted clusters.
  #'
  #' @export
  #'
  records <- initialize_records(dataset_init, init_population)
  results <- feve_main(init_population, dataset_init, records, params, figures, sheets)
  return(results)
}

resume_feve <- function(dataset_init, records, params, figures=TRUE, sheets=TRUE) {
  #' Conduct a clustering analysis with the fEVE framework, starting from the last pending population.
  #'
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #' @param sheets a boolean that indicates if the results of the clustering recursion should be saved in Excel sheets.
  #'
  #' @return a named list, with two elements: `records` and `preds`.
  #' `records` is a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' `preds` is a named vector associating samples to their predicted clusters.
  #'
  #' @export
  #'
  population <- get_pending_population(records$meta)
  results <- feve_main(population, dataset_init, records, params, figures, sheets)
  return(results)
}
