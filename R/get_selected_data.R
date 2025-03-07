"Functions called to set-up the data of a fEVE clustering recursion.

	2025/03/05 @yanisaspic"

get_selected_dataset <- function(dataset, params) {
  #' Get a dataset with selected features.
  #'
  #' @param dataset a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a dataset, with selected features.
  #'
  selected_features <- params$selected_features_strategy(dataset, params)
  selected_dataset <- dataset[selected_features,]
  return(selected_dataset)
}

get_selected_SeuratObject <- function(selected_dataset, params) {
  #' Get a SeuratObject from a dataset with selected features.
  #'
  #' @param selected_dataset a dataset, with selected features.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a SeuratObject, on which the function ScaleData() of Seurat has been applied already.
  #'
  #' @import Seurat
  #'
  SeurObj <- Seurat::CreateSeuratObject(selected_dataset)
  Seurat::VariableFeatures(SeurObj) <- rownames(selected_dataset)
  SeurObj <- Seurat::NormalizeData(SeurObj)
  SeurObj <- Seurat::ScaleData(SeurObj, features=Seurat::VariableFeatures(SeurObj))
  SeurObj <- Seurat::RunPCA(SeurObj, features=Seurat::VariableFeatures(SeurObj), seed.use=params$random_state)
  SeurObj <- Seurat::RunUMAP(SeurObj, features=Seurat::VariableFeatures(SeurObj), seed.use=params$random_state)
  return(SeurObj)
}

get_selected_data <- function(population, dataset_init, SeuratObject_init, records, params, figures) {
  #' Extract a data subset with a specific population, and its most variable features.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param SeuratObject_init a SeuratObject generated from dataset_init, on which
  #' the function RunUMAP() of Seurat has been applied already.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #'
  #' @return a named list, with two names: `dataset` and `SeuratObject`.
  #'
  #' @import Seurat
  #' @import glue
  #' @import grDevices
  #' @import qpdf
  #'
  #' @export
  #'
  samples_of_population <- get_samples_of_population(population, records$samples)
  if (length(samples_of_population) >= params$minimum_samples) {
    tmp <- dataset_init[, samples_of_population]
    selected_dataset <- get_selected_dataset(tmp, params)
    selected_SeuratObject <- get_selected_SeuratObject(selected_dataset, params)
    selected_data <- list(dataset=selected_dataset,
                          SeuratObject=selected_SeuratObject)}
  else {selected_data <- list()}

  if (figures) {
    plot <- draw_selected_data(population, SeuratObject_init, records)
    grDevices::pdf(file = glue::glue("{params$figures_path}/{population}_selected_data.pdf"))
    print(plot)
    grDevices::dev.off()
  }
  return(selected_data)
}

draw_selected_data <- function(population, SeuratObject_init, records) {
  #' Get a U-MAP plot representing the pool of samples used in the clustering recursion.
  #'
  #' @param population a character. It corresponds to the cell population that scEVE will attempt to cluster.
  #' @param SeuratObject_init a SeuratObject generated from dataset_init, on which
  #' the function RunUMAP() of Seurat has been applied already.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @return a plot.
  #'
  #' @import viridis
  #' @import forcats
  #' @import assertthat
  #' @import ggplotify
  #' @import SCpubr
  #' @import ggplot2
  #'
  #' @export
  #'
  samples_of_population <- get_samples_of_population(population, records$samples)
  plot <- SCpubr::do_DimPlot(SeuratObject_init, cells.highlight=samples_of_population) +
    ggplot2::ggtitle(population) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.background=ggplot2::element_rect(fill="lightgrey"),
                   axis.title=ggplot2::element_blank(), legend.position="bottom")
  return(plot)
}
