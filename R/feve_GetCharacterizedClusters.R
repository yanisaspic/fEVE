"Functions called to identify characterized clusters, i.e. clusters with sufficient marker features.

	2025/04/28 @yanisaspic"

feve_GetCharacterizedClusters <- function(clusters, selected_data, params, features_threshold=10) {
  #' Filter out uncharacterized clusters, i.e. clusters with too little marker features.
  #'
  #' @param clusters a list where every element is a pool of samples.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features` and `robustness`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param features_threshold the minimal number of marker features expected in a characterized cluster.
  #'
  #' @return a list where every element is a pool of samples.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features` and `robustness`.
  #'
  #' @export
  #'
  is_characterized <- function(cluster) {length(cluster$features) >= features_threshold}
  characterized_clusters <- Filter(f=is_characterized, x=clusters)
  return(characterized_clusters)
}