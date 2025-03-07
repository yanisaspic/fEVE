"Functions called to report the results of a fEVE clustering recursion.

	2025/03/05 @yanisaspic"

get_drawings_paths <- function(population, params) {
  #' Get the paths leading to each informative figure drawn during the clustering recursion.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a vector of paths.
  #' 
  #' @import glue
  #'
  drawings_paths <- c()
  for (drawing in c("selected_data", "base_clusters", "robust_clusters", "characteristic_features")) {
    filename <- glue::glue("{params$figures_path}/{population}_{drawing}.pdf")
    if (file.exists(filename)) {drawings_paths <- c(drawings_paths, filename)}}
  return(drawings_paths)
}

merge_drawings <- function(population, params) {
  #' Merge the informative figures drawn during the clustering recursion in a single .pdf file.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @import glue
  #' @import qpdf
  #'
  drawings_paths <- get_drawings_paths(population, params)
  qpdf::pdf_combine(input=drawings_paths, output=glue::glue("{params$figures_path}/{population}.pdf"))
  unlink(drawings_paths)
}

report_samples <- function(characterized_clusters, selected_data, records, params) {
  #' Report the composition of the characterized clusters in the corresponding records sheet.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of samples.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features`, `robustness` and `specific_features`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are characterized clusters, and cluster memberships are reported in the table.
  #'
  cluster_memberships <- params$cluster_memberships_strategy(characterized_clusters, selected_data, params)
  records$samples <- merge(records$samples, cluster_memberships, by="row.names", all=TRUE)
  records$samples <- transform(records$samples, row.names=Row.names, Row.names=NULL)
  records$samples[is.na(records$samples)] <- 0
  return(records$samples)
}

report_features <- function(characterized_clusters, records) {
  #' Report the characteristic features of the characterized clusters in the corresponding records sheet.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of samples.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features`, `robustness` and `specific_features`.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @return a data.frame associating predicted populations to their characteristic features.
  #' Its rows are genes, its columns are predicted populations,
  #' and the strength of the characterization (e.g. log2-transformed pvalues) are reported in the table.
  #'
  #' @import stats
  #'
  get_column <- function(cluster) {
    column <- stats::setNames(data.frame(cluster$features), cluster$label)
    column[setdiff(rownames(records$features), rownames(column)),] <- 0
    return(column)}
  columns <- lapply(X=characterized_clusters, FUN=get_column)
  for (c in columns) {
    records$features <- merge(records$features, c, by="row.names", all=TRUE)
    records$features <- transform(records$features, row.names=Row.names, Row.names=NULL)}
  return(records$features)
}

report_methods <-function(characterized_clusters, records) {
  #' Report the clustering methods used to predict the characterized clusters in the corresponding records sheet.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of samples.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features`, `robustness` and `specific_features`.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @return a data.frame associating predicted populations to the clustering methods leveraged to predict them.
  #' Its rows are clustering methods, its columns are predicted populations, and binary values are reported in the table.
  #'
  methods <- sapply(X=characterized_clusters, FUN="[[", "clustering_methods")
  data <- utils::stack(methods) %>% dplyr::rename(methods=values, clusters=ind)
  methods.recursion <- table(data$methods, data$clusters)
  methods.recursion <- as.data.frame.matrix(methods.recursion)
  records$methods <- merge(records$methods, methods.recursion, by="row.names", all=TRUE)
  records$methods <- transform(records$methods, row.names=Row.names, Row.names=NULL)
  records$methods[is.na(records$methods)] <- 0
  return(records$methods)
}

report_metadata <- function(population, characterized_clusters, records) {
  #' Report the metadata regarding characterized clusters in the corresponding records sheet.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param characterized_clusters a list where every element is a characterized pool of samples.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features`, `robustness` and `specific_features`.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @return a data.frame associating predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #'
  get_metadata.cluster <- function(cluster) {
    samples_of_population <- get_samples_of_population(cluster$label, records$samples)
    metadata.cluster <- c(size=length(samples_of_population), robustness=cluster$robustness,
                          parent=population, clustering_status="PENDING")
    return(metadata.cluster)}
  metaselected_data <- lapply(X=characterized_clusters, FUN=get_metadata.cluster)
  metaselected_data <- do.call(rbind, metaselected_data)
  records$meta <- rbind(records$meta, metaselected_data)
  for (col in c("robustness", "size")) {records$meta[, col] <- as.numeric(records$meta[, col])}
  return(records$meta)
}

get_recorded_recursion <- function(population, characterized_clusters, selected_data, records, params) {
  #' Report multiple information regarding the characterized clusters predicted during a clustering recursion.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param characterized_clusters a list where every element is a characterized pool of samples.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features`, `robustness` and `specific_features`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param records a named list, with three data.frames: `samples`, `features` and `meta`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @export
  #'
  records$samples <- report_samples(characterized_clusters, selected_data, records, params)
  records$features <- report_features(characterized_clusters, records)
  records$meta <- report_metadata(population, characterized_clusters, records)
  records$methods <- report_methods(characterized_clusters, records)
  return(records)
}
