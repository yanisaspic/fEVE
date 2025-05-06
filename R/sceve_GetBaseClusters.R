"Functions used to get base clusters.

	2025/03/06 @yanisaspic"

use_Seurat <- function(selected_SeuratObject, params) {
  #' Predict clusters with the Seurat method.
  #'
  #' @param selected_SeuratObject a SeuratObject, on which the function ScaleData()
  #' of Seurat has been applied already.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import Seurat
  #'
  if (!"umap" %in% names(selected_SeuratObject@reductions)) {
    selected_SeuratObject <- Seurat::RunPCA(selected_SeuratObject,
                                            features=Seurat::VariableFeatures(selected_SeuratObject),
                                            seed.use=params$random_state)}
  selected_SeuratObject <- Seurat::FindNeighbors(selected_SeuratObject,
                                                 features=Seurat::VariableFeatures(selected_SeuratObject))
  selected_SeuratObject <- Seurat::FindClusters(selected_SeuratObject,
                                                random.seed=params$random_state)
  predictions <- Seurat::Idents(selected_SeuratObject)
  return(predictions)
}

use_monocle3 <- function(selected_SeuratObject, params) {
  #' Predict clusters with the monocle3 method.
  #'
  #' @param selected_SeuratObject a SeuratObject, on which the function ScaleData()
  #' of Seurat has been applied already.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import monocle3
  #' @import R.utils
  #' @import Seurat
  #' @import SeuratWrappers
  #'
  if (!"umap" %in% names(selected_SeuratObject@reductions)) {
    selected_SeuratObject <- Seurat::RunUMAP(selected_SeuratObject,
                                             features=Seurat::VariableFeatures(selected_SeuratObject),
                                             seed.use=params$random_state)}
  CDSObject <- SeuratWrappers::as.cell_data_set(selected_SeuratObject)
  CDSObject <- monocle3::cluster_cells(CDSObject, random_seed=params$random_state)
  predictions <- CDSObject@clusters@listData$UMAP$clusters
  return(predictions)
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

use_SHARP <- function(selected_dataset, params) {
  #' Predict clusters with the SHARP method.
  #'
  #' @param selected_dataset a scRNA-seq dataset of raw count expression, with selected genes.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import clues
  #' @import SHARP
  #'
  results <- SHARP::SHARP(scExp=selected_dataset, exp.type="count",
                          n.cores = 1, rN.seed=params$random_state)
  predictions <- get_formatted_predictions(samples=colnames(selected_dataset),
                                           predictions=results$pred_clusters)
  return(predictions)
}

use_densityCut <- function(logtpm_dataset, params) {
  #' Predict clusters with the densityCut method.
  #'
  #' @param logtpm_dataset a scRNA-seq dataset of log-transformed transcripts per millions (tpm),
  #' with selected genes.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  #' @import densitycut
  #'
  set.seed(params$random_state)
  data <- t(logtpm_dataset) # densityCut expects cells as rows and genes as columns.
  results <- densitycut::DensityCut(t(logtpm_dataset), show.plot = FALSE)
  predictions <- get_formatted_predictions(samples=colnames(logtpm_dataset), predictions=results$cluster)
  return(predictions)
}

get_data_input <- function(method, data) {
  #' Get the data input relevant to a clustering method.
  #'
  #' @param method a clustering method used by scEVE.
  #' Valid methods include: densityCut, monocle3, Seurat and SHARP.
  #' @param data a named list, with three names: `dataset`, `SeuratObject` and `logtpm`.
  #' They correspond to the scRNA-seq expression matrix of a specific cell population, its SeuratObject and its log2-transformed TPM values.
  #'
  #' @return one of `data$dataset`, `data$SeuratObject` or `data$logtpm`.
  #'
  data_inputs <- list(Seurat=data$SeuratObject,
                      monocle3=data$SeuratObject,
                      SHARP=data$dataset,
                      densityCut=data$logtpm)
  input <- data_inputs[[method]]
  return(input)
}

sceve_GetBaseClusters <- function(selected_data, params,
                                  clustering_methods=c("densityCut", "monocle3", "Seurat", "SHARP")) {
  #' Apply multiple clustering methods on scRNA-seq data to predict base clusters.
  #'
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param clustering_methods a vector of valid clustering methods.
  #'
  #' @return a data.frame associating cells to their predicted clusters.
  #' Its rows are cells, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @import glue
  #' @import scater
  #'
  #' @export
  #'
  logtpm_dataset <- log2(scater::calculateTPM(selected_data$dataset) + 1)
  selected_data[["logtpm"]] <- logtpm_dataset

  get_predictions_method <- function(method) {
    f <- get(glue::glue("use_{method}"))
    predictions <- f(get_data_input(method, selected_data), params)
    return(predictions)
  }
  base_clusters <- sapply(X=clustering_methods, FUN=get_predictions_method)
  gc()

  rename_clusters_method <- function(method) {
    rename_prediction <- function(prediction) {glue::glue("{method}_{prediction}")}
    column <- sapply(X=base_clusters[, method], FUN=rename_prediction)
    return(column)
  }
  base_clusters <- sapply(X=clustering_methods, FUN=rename_clusters_method)
  return(base_clusters)
}
