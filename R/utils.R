"Some useful functions developed to interact with fEVE data.

	2025/05/09 @yanisaspic"

initialize_records <- function(dataset_init, init_population) {
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
  #' @param init_population a character (without `.`).
  #'
  #' @return a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @export
  #'
  samples <- data.frame(as.numeric(rep(1, ncol(dataset_init))), row.names=colnames(dataset_init))
  features <- data.frame(as.numeric(rep(0, nrow(dataset_init))), row.names=rownames(dataset_init))
  meta <- data.frame(size=as.numeric(ncol(dataset_init)), robustness=0, parent=NA,
                     clustering_status="PENDING", row.names=init_population)
  methods <- data.frame()
  records <- list(samples=samples, features=features, meta=meta, methods=methods)
  for (sheet in c("samples", "features")) {colnames(records[[sheet]]) <- init_population}
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

get_pending_population <- function(records_meta) {
  #' Get a population for which no fEVE clustering recursion has been attempted.
  #'
  #' @param records_meta a data.frame associating predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #'
  #' @return a character.
  #'
  #' @export
  #'
  pending_populations <- rownames(records_meta[records_meta$clustering_status=="PENDING", ])
  population <- pending_populations[1]
  return(population)
}

get_resolution <- function(population) {
  #' Get the resolution level of a predicted population.
  #'
  #' @param population a character. It corresponds to a population predicted by fEVE (e.g. C.2).
  #'
  #' @return an integer.
  #'
  resolution <- length(strsplit(population, split=".", fixed=TRUE)[[1]])
  return(resolution)
}

get_populations_at_resolution <- function(resolution, records_samples) {
  #' Get all the populations at a specific resolution.
  #'
  #' The root population is resolution 1, and its children populations are resolution 2, etc.
  #'
  #' @param resolution an integer.
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a vector of population labels.
  #'
  populations <- colnames(records_samples)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  populations_at_resolution <- populations[resolutions==resolution]
  return(populations_at_resolution)
}

get_maximum_resolution <- function(records_samples) {
  #' Get the maximum clustering resolution attained in a fEVE clustering analysis.
  #'
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are populations. The cell values range from 0 to 1.
  #'
  #' @return an integer.
  #'
  populations <- colnames(records_samples)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  maximum_resolution <- max(resolutions)
  return(maximum_resolution)
}

get_samples_of_population <- function(population, records_samples) {
  #' Get every sample belonging in a predicted population.
  #'
  #' @param population a character. It corresponds to a population predicted by fEVE (e.g. C.2).
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a vector of samples.
  #'
  init_population <- colnames(records_samples)[1]
  if (population==init_population) {return(rownames(records_samples))}
  # all samples belong to the root population

  cell_is_in_population <- function(cell_membership) {
    max_membership <- max(cell_membership)
    main_population <- names(which.max(cell_membership))
    is_in_population <- (main_population == population) &
      (max_membership > 0)
  }

  resolution <- get_resolution(population)
  records_samples <- records_samples[, get_populations_at_resolution(resolution, records_samples)]
  samples_are_in_population <- apply(X=records_samples, MARGIN=1, FUN=cell_is_in_population)
  samples_of_population <- rownames(records_samples[samples_are_in_population, ])
  return(samples_of_population)
}

get_records <- function(sheets_path) {
  #' Load the records of a fEVE clustering analysis.
  #'
  #' @param sheets_path a path where result sheets are stored.
  #'
  #' @return a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @import openxlsx
  #'
  #' @export
  #'
  sheets_names <- openxlsx::getSheetNames(sheets_path)
  get_sheet <- function(sheets_name) {openxlsx::read.xlsx(sheets_path, sheet=sheets_name, rowNames=TRUE)}
  records <- sapply(X=sheets_names, FUN=get_sheet)
  return(records)
}

get_leaf_clusters <- function(records_samples) {
  #' Get a named vector associating each sample to its most informative cluster label.
  #'
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a named vector associating samples to their most informative cluster labels.
  #'
  #' @import stats
  #'
  #' @export
  #'
  if (ncol(records_samples)==1) {
    init_population <- colnames(records_samples)[1]
    leaf_clusters <- stats::setNames(object=rep(init_population, nrow(records_samples)), nm=rownames(records_samples))
    leaf_clusters <- leaf_clusters[order(names(leaf_clusters))]
    return(leaf_clusters)}

  get_labels_at_resolution <- function(resolution) {
    populations_at_resolution <- get_populations_at_resolution(resolution, records_samples)
    get_labels.population <- function(population) {
      samples_of_population <- get_samples_of_population(population, records_samples)
      labels.population <- stats::setNames(rep(population, length(samples_of_population)), samples_of_population)
      return(labels.population)}
    labels_at_resolution <- sapply(X=populations_at_resolution, FUN=get_labels.population)
    labels_at_resolution <- unlist(unname(labels_at_resolution))
    return(labels_at_resolution)
  }

  # get the labels of every sample with a bottom-up approach, where sample labels are successively added
  # from the maximum resolution to the minimal one (i.e. the label 'C').
  labels <- lapply(X=get_maximum_resolution(records_samples):2, FUN=get_labels_at_resolution)
  labels <- unlist(labels)
  labels <- labels[!duplicated(names(labels))]
  labels <- labels[order(names(labels))]
  return(labels)
}

get_leaf_clusters_at_resolution <- function(resolution, records_samples) {
  #' Get a named vector associating each sample to its most informative cluster label, at a given maximum resolution.
  #'
  #' @param resolution an integer.
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a named vector associating samples to their most informative cluster labels, at a given maximum resolution.
  #'
  f <- function(resolution) {get_populations_at_resolution(resolution, records_samples)}
  tmp <- sapply(X=1:resolution, FUN=f)
  if (length(tmp) > 1) {tmp <- do.call(c, tmp)}
  data_at_resolution <- records_samples[, tmp, drop=FALSE]
  leaf_clusters_at_resolution <- get_leaf_clusters(data_at_resolution)
  return(leaf_clusters_at_resolution)
}