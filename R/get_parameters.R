"Function called to load a subset of parameters

	2025/05/09 @yanisaspic"

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

  # default fEVE instance for multiple -omics data ______________________________
  params[["fEVE"]] <- list(random_state=1, minimum_samples=51, # npcs=50 for Seurat
                          figures_path="./fEVE", sheets_path="./fEVE/records.xlsx",
                          selected_features_strategy=feve_GetSelectedFeatures,
                          base_clusters_strategy=feve_GetBaseClusters,
                          characteristic_features_strategy=feve_GetCharacteristicFeatures,
                          characterized_clusters_strategy=feve_GetCharacterizedClusters,
                          cluster_memberships_strategy=feve_HardClustering)

  # scEVE instance for single-cell transcriptomics data __________________________
  params[["scEVE"]] <- list(random_state=1, minimum_samples=100, # sncells=100 for SHARP
                          figures_path="./scEVE", sheets_path="./scEVE/records.xlsx",
                          selected_features_strategy=sceve_GetSelectedFeatures,
                          base_clusters_strategy=sceve_GetBaseClusters,
                          characteristic_features_strategy=sceve_GetCharacteristicFeatures,
                          characterized_clusters_strategy=sceve_GetCharacterizedClusters,
                          cluster_memberships_strategy=feve_HardClustering)

  return(params[[params_label]])
}