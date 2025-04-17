"Functions used to get base clusters.

	2025/03/18 @yanisaspic"

get_transformed_dataset <- function(filtered_dataset) {
  #' Get a log-transformed and normalized dataset.
  #'
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #'
  #' @return a log-transformed and normalized -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #'
  #' @import SNFtool
  #'
  dataset <- log(filtered_dataset + 1)
  dataset <- SNFtool::standardNormalization(t(dataset)) # normalize features in columns
  return(t(dataset))
}

get_formatted_predictions <- function(samples, predictions) {
  #' Format the cluster predictions of a method to obtain a standard output.
  #'
  #' @param samples a vector of sample names.
  #' @param predictions a vector of cluster predictions.
  #'
  #' @return a named factor that associates each sample to their cluster prediction.
  #'
  #' @import stats
  #'
  predictions <- stats::setNames(predictions, samples)
  predictions <- factor(predictions)
  return(predictions)
}

use_NEMO <- function(filtered_dataset, params) {
  #' Predict clusters with the NEMO method.
  #'
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each patient to their cluster prediction.
  #'
  #' @import NEMO
  #'
  set.seed(params$random_state)
  clusters <- NEMO::nemo.clustering(list(filtered_dataset))
  predictions <- get_formatted_predictions(samples=names(clusters), predictions=unname(clusters))
  return(predictions)
}

use_PINS <- function(filtered_dataset, params) {
  #' Predict clusters with the PINS method.
  #'
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each patient to their cluster prediction.
  #'
  #' @import PINSPlus
  #'
  set.seed(params$random_state)
  output <- PINSPlus::PerturbationClustering(data=t(filtered_dataset), ncore=4, kMax=15)
  predictions <- get_formatted_predictions(samples=names(output$cluster), predictions=unname(output$cluster))
  return(predictions)
}

use_SNF <- function(filtered_dataset, params) {
  #' Predict clusters with the SNF method.
  #'
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each patient to their cluster prediction.
  #'
  #' @import SNFtool
  #'
  set.seed(params$random_state)
  distances <- SNFtool::dist2(as.matrix(t(filtered_dataset)), as.matrix(t(filtered_dataset)))
  sigma <- 0.5
  k <- round(1/10 * ncol(filtered_dataset))  # Number of nearest neighbors = 1/10 of the nb of samples
  aff <- SNFtool::affinityMatrix(distances, k, sigma)
  num_clusters <- SNFtool::estimateNumberOfClustersGivenGraph(aff, 2:15)[[3]]  # between 2 and 15 clusters, rotation cost method
  clusters <- SNFtool::spectralClustering(aff, num_clusters)
  predictions <- get_formatted_predictions(samples=colnames(filtered_dataset), predictions=clusters)
  return(predictions)
}

use_kMeans <- function(filtered_dataset, params) {
  #' Predict clusters with the SNF method.
  #'
  #' @param filtered_dataset an -omics dataset, with variant features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each patient to their cluster prediction.
  #'
  #' @import cluster
  #' @import stats
  #'
  set.seed(params$random_state)
  sil <- rep(0, 20)
  clusterings <- list()

  # repeat k-means for k 1:20 and extract silhouette:
  for (i in 2:20) {
    k1to20 <- stats::kmeans(t(filtered_dataset), centers = i)
    clusterings[[i]] <- k1to20
    ss <- cluster::silhouette(k1to20$cluster, dist(t(filtered_dataset)))
    sil[i] <- mean(ss[, 3])}

  kopt = which.max(sil)
  clusters = clusterings[[kopt]]$cluster
  predictions <- get_formatted_predictions(samples=names(clusters), predictions=unname(clusters))
  return(predictions)
}

breve_GetBaseClusters <- function(selected_data, params,
                                    clustering_methods=c("NEMO", "PINS", "SNF", "kMeans")) {
  #' Apply multiple clustering methods on an -omics dataset to predict base clusters.
  #'
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param clustering_methods a vector of valid clustering methods.
  #'
  #' @return a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @import glue
  #'
  #' @export
  #'
  selected_data$dataset <- get_transformed_dataset(selected_data$dataset)
  get_predictions_method <- function(method) {
    f <- get(glue::glue("use_{method}"))
    predictions <- f(selected_data$dataset, params)
    return(predictions)}
  base_clusters <- sapply(X=clustering_methods, FUN=get_predictions_method)
  gc()

  rename_clusters_method <- function(method) {
    rename_prediction <- function(prediction) {glue::glue("{method}_{prediction}")}
    column <- sapply(X=base_clusters[, method], FUN=rename_prediction)
    return(column)}

  base_clusters <- sapply(X=clustering_methods, FUN=rename_clusters_method)
  return(base_clusters)
}
