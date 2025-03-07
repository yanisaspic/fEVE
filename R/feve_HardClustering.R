"Functions used to report the results of a clustering recursion as hard clusters.

	2025/03/06 @yanisaspic"

feve_HardClustering <- function(characterized_clusters, selected_data, params) {
  #' Get a table reporting the composition of the characterized clusters.
  #'
  #' By calling this function, the samples are hard clustered,
  #' i.e. their membership to their assigend cluster is 1, and their membership to the other clusters is 0.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of samples.
  #' The elements are named lists, with seven names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features`, `robustness` and `specific_features`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are characterized clusters, and cluster memberships are reported in the table.
  #'
  #' @import dplyr
  #' @importFrom rlang .data
  #' @import utils
  #'
  #' @export
  #'
  samples <- sapply(X=characterized_clusters, FUN="[[", "samples")
  data <- utils::stack(samples) %>% dplyr::rename(samples=values, clusters=ind)
  samples_recursion <- table(data$samples, data$clusters)
  samples_recursion <- as.data.frame.matrix(samples_recursion)
  return(samples_recursion)
}