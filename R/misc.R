"Miscellaneous functions called multiple times in the fEVE framework.

	2025/03/05 @yanisaspic"

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

get_populations_at_resolution <- function(resolution, records.samples) {
  #' Get all the populations at a specific resolution.
  #'
  #' The root population is resolution 1, and its children populations are resolution 2, etc.
  #'
  #' @param resolution an integer.
  #' @param records.samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a vector of population labels.
  #'
  populations <- colnames(records.samples)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  populations_at_resolution <- populations[resolutions==resolution]
  return(populations_at_resolution)
}

get_maximum_resolution <- function(records.samples) {
  #' Get the maximum clustering resolution attained in a fEVE clustering analysis.
  #'
  #' @param records.samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and and its columns are populations. The cell values range from 0 to 1.
  #'
  #' @return an integer.
  #'
  populations <- colnames(records.samples)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  maximum_resolution <- max(resolutions)
  return(maximum_resolution)
}

get_samples_of_population <- function(population, records.samples) {
  #' Get every sample belonging in a predicted population.
  #'
  #' @param population a character. It corresponds to a population predicted by fEVE (e.g. C.2).
  #' @param records.samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a vector of samples.
  #'
  #' @export
  #'
  init_population <- colnames(records.samples)[1]
  if (population==init_population) {return(rownames(records.samples))}
  # all samples belong to the root population

  cell_is_in_population <- function(cell_membership) {
    max_membership <- max(cell_membership)
    main_population <- names(which.max(cell_membership))
    is_in_population <- (main_population == population) &
      (max_membership > 0)
  }

  resolution <- get_resolution(population)
  records.samples <- records.samples[, get_populations_at_resolution(resolution, records.samples)]
  samples_are_in_population <- apply(X=records.samples, MARGIN=1, FUN=cell_is_in_population)
  samples_of_population <- rownames(records.samples[samples_are_in_population, ])
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

get_leaf_clusters <- function(records.samples) {
  #' Get a named vector associating each sample to its most informative cluster label.
  #'
  #' @param records.samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a named vector associating samples to their most informative cluster labels.
  #'
  #' @import stats
  #'
  #' @export
  #'
  if (ncol(records.samples)==1) {
    init_population <- colnames(records.samples)[1]
    leaf_clusters <- stats::setNames(object=rep(init_population, nrow(records.samples)), nm=rownames(records.samples))
    leaf_clusters <- leaf_clusters[order(names(leaf_clusters))]
    return(leaf_clusters)}

  get_labels.resolution <- function(resolution) {
    populations_at_resolution <- get_populations_at_resolution(resolution, records.samples)
    get_labels.population <- function(population) {
      samples_of_population <- get_samples_of_population(population, records.samples)
      labels.population <- stats::setNames(rep(population, length(samples_of_population)), samples_of_population)
      return(labels.population)}
    labels.resolution <- sapply(X=populations_at_resolution, FUN=get_labels.population)
    labels.resolution <- unlist(unname(labels.resolution))
    return(labels.resolution)
  }

  # get the labels of every sample with a bottom-up approach, where sample labels are successively added
  # from the maximum resolution to the minimal one (i.e. the label 'C').
  labels <- lapply(X=get_maximum_resolution(records.samples):2, FUN=get_labels.resolution)
  labels <- unlist(labels)
  labels <- labels[!duplicated(names(labels))]
  labels <- labels[order(names(labels))]
  return(labels)
}

get_leaf_clusters.resolution <- function(resolution, records.samples) {
  #' Get a named vector associating each sample to its most informative cluster label, at a given maximum resolution.
  #'
  #' @param resolution an integer.
  #' @param records.samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and and its columns are population. The cell values range from 0 to 1.
  #'
  #' @return a named vector associating samples to their most informative cluster labels, at a given maximum resolution.
  #'
  data.resolution <- records.samples[, get_populations_at_resolution(resolution, records.samples)]
  leaf_clusters.resolution <- get_leaf_clusters(data.resolution)
  return(leaf_clusters.resolution)
}
