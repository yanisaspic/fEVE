"Main functions called to condut a fEVE clustering analysis.

	2025/03/06 @yanisaspic"

get_parameters <- function(params_label) {
  #' Get parameters to conduct a fEVE clustering analysis.
  #'
  #' @param params_label a character.
  #'
  #' These parameters are:
  #' `random_state` the integer seed used to have deterministic results.
  #' `minimum_samples` the minimum number of samples expected in a population before attempting to sub-cluster it.
  #' `figures_path` & `sheets_path` the default paths where figures and result sheets are stored, respectively.
  #' `selected_features_strategy` a function called to select a limited pool of features for a clustering recursion.
  #' `base_clusters_strategy` a function called to predict base clusters with multiple clustering methods.
  #' `characteristic_features_strategy` a function called to predict the characteristic features of every meta-cluster.
  #' `characterized_clusters_strategy` a function called to identify characterized clusters.
  #' `cluster_memberships_strategy` a function called to assign ambiguous samples to their clusters (i.e. hard or soft-clustering).
  #'
  #' Detailed information regarding the strategy parameters are available in the vignette of the package.
  #'
  #' @return a list of parameters.
  #'
  #' @export
  #'
  params <- list()

  # scEVE for single-cell transcriptomics data _________________________________
  params[["scEVE"]] <- list(random_state=1, minimum_samples=100, # sncells=100 for SHARP
                          figures_path="./scEVE", sheets_path="./scEVE/records.xlsx",
                          selected_features_strategy=sceve_GetSelectedFeatures,
                          base_clusters_strategy=sceve_GetBaseClusters,
                          characteristic_features_strategy=sceve_GetCharacteristicFeatures,
                          characterized_clusters_strategy=sceve_GetCharacterizedClusters,
                          cluster_memberships_strategy=feve_HardClustering)

  params[["brEVE"]] <- list(random_state=1, minimum_samples=51, # npcs=50 for Seurat
                          figures_path="./brEVE", sheets_path="./breve/records.xlsx",
                          selected_features_strategy=breve_GetSelectedFeatures,
                          base_clusters_strategy=breve_GetBaseClusters,
                          characteristic_features_strategy=breve_GetCharacteristicFeatures,
                          characterized_clusters_strategy=sceve_GetCharacterizedClusters,
                          cluster_memberships_strategy=feve_HardClustering)

  return(params[[params_label]])
}

initialize_records <- function(dataset_init) {
  #' Get a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' - `samples` associates samples to their predicted populations.
  #' Its rows are samples, its columns are predicted populations, and cluster memberships are reported in the table.
  #' - `features` associates predicted populations to their characteristic features.
  #' Its rows are features, its columns are predicted populations, and characterization powers are reported in the table.
  #' - `meta` associates predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #' - `methods` associates predicted populations to the clustering methods leveraged to predict them.
  #' Its rows are clustering methods, its columns are predicted populations, and binary values are reported in the table.
  #'
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #'
  #' @return a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @export
  #'
  samples <- data.frame(C=as.numeric(rep(1, ncol(dataset_init))), row.names=colnames(dataset_init))
  features <- data.frame(C=as.numeric(rep(0, nrow(dataset_init))), row.names=rownames(dataset_init))
  meta <- data.frame(size=as.numeric(ncol(dataset_init)), robustness=0, parent=NA,
                     clustering_status="PENDING", row.names="C")
  methods <- data.frame()
  records <- list(samples=samples, features=features, meta=meta, methods=methods)
  return(records)
}

get_SeuratObject_init <- function(dataset_init) {
  #' Get a SeuratObject from a dataset, without selected features.
  #'
  #' This function is used once prior to a fEVE clustering analysis in order to draw the
  #' extracted data at each clustering recursion (cf. `feve::draw_extracted_data()`).
  #'
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #'
  #' @return a SeuratObject, on which the function RunUMAP() of Seurat has been applied already.
  #'
  #' @import Seurat
  #'
  #' @export
  #'
  SeuratObject_init <- Seurat::CreateSeuratObject(dataset_init)
  SeuratObject_init <- Seurat::FindVariableFeatures(SeuratObject_init)
  SeuratObject_init <- Seurat::NormalizeData(SeuratObject_init)
  SeuratObject_init <- Seurat::ScaleData(SeuratObject_init,
                                         features=Seurat::VariableFeatures(SeuratObject_init))
  SeuratObject_init <- Seurat::RunUMAP(SeuratObject_init,
                                       features=Seurat::VariableFeatures(SeuratObject_init),
                                       seed.use=1)
  return(SeuratObject_init)
}

get_pending_population <- function(records) {
  #' Get a population for which no fEVE clustering recursion has been attempted.
  #'
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @return a character.
  #'
  #' @export
  #'
  meta <- records$meta
  pending_populations <- rownames(meta[meta$clustering_status=="PENDING", ])
  population <- pending_populations[1]
  return(population)
}

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
  #' @import openxlsx
  #'
  if (figures) {
    dir.create(params$figures_path)
    SeuratObject_init <- get_SeuratObject_init(dataset_init)}
  else {SeuratObject_init <- NA}

  while (!is.na(population)) {
    print(population)
    records <- feve_recursion(population, dataset_init, SeuratObject_init, records, params, figures, sheets)
    population <- get_pending_population(records)}

  feature_is_characteristic <- function(feature) {sum(feature) > 0}
  records$features <- records$features[apply(X=records$features, MARGIN=1, FUN=feature_is_characteristic),]
  if (sheets) {openxlsx::write.xlsx(records, params$sheets_path, rowNames=TRUE)}

  results <- list(records=records, preds=factor(get_leaf_clusters(records$samples)))
  return(results)
}

feve <- function(dataset_init, params, figures=TRUE, sheets=TRUE) {
  #' Conduct a clustering analysis with the fEVE framework, starting from the initial population C.
  #'
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #' @param sheets a boolean that indicates if the results of the clustering recursion should be saved in Excel sheets.
  #'
  #' @return a named list, with two elements: `records` and `preds`.
  #' `records` is a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' `preds` is a named factor associating samples to their predicted clusters.
  #'
  #' @export
  #'
  records <- initialize_records(dataset_init)
  population <- "C"
  results <- feve_main(population, dataset_init, records, params, figures, sheets)
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
  #' `preds` is a named factor associating samples to their predicted clusters.
  #'
  #' @export
  #'
  population <- get_pending_population(records)
  results <- feve_main(population, dataset_init, records, params, figures, sheets)
  return(results)
}
