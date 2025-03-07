"Function called to load scRNA-seq datasets.

	2025/03/06 @yanisaspic"

sceve_GetDatasets <- function() {
  #' Get information regarding the datasets available on TMExplorer for scEVE.
  #'
  #' @return a data.frame associating each dataset to its sequencing protocol,
  #' its number of cells, clusters and genes, as well as its year of publication,
  #' its accession number and its associated doi.
  #'
  #' @export
  #'
  characteristics <- c("sequencing_protocol", "n_cells", "n_clusters", "n_genes",
                       "year", "accession_number", "doi")
  metadata <- c("SMARTer (Fluidigm C1)", 364, 7, 57241, 2017, "GSE81861", "10.1038/ng.3818",
                "SMART-Seq2", 3589, 7, 23460, 2017, "GSE84465", "10.1016/j.celrep.2017.10.030",
                "SMART-Seq2", 6879, 9, 23686, 2018, "GSE115978", "10.1016/j.cell.2018.09.006",
                "10x Genomics", 18456, 18, 23580, 2020, "GSE125969", "10.1016/j.celrep.2020.108023",
                "Seq-Well", 22600, 17, 27899, 2018, "GSE116256", "10.1016/j.cell.2019.01.031",
                "10x Genomics", 51775, 17, 22180, 2018, "E-MTAB-6149,   E-MTAB-6653,", "10.1038/s41591-018-0096-5",
                "10x Genomics", 57530, 10, 24005, 2019, "CRA001160", "10.1038/s41422-019-0195-y")
  datasets <- c("Li_HumCRC_b", "Darmanis_HumGBM", "JerbyArnon_HumMLM", "Gillen_HumEPN",
                "VanGalen_HumAML", "Lambrechts_HumNSCLC", "Peng_HumPDAC")
  output <- matrix(metadata, nrow=length(datasets), ncol=length(characteristics), byrow=TRUE)
  output <- as.data.frame.matrix(output, row.names=datasets)
  colnames(output) <- characteristics
  return(output)
}

get_cell_ids <- function(labels) {
  #' Get a vector of unique cell ids according to the label of each cell.
  #'
  #' @param labels a vector of cell types.
  #'
  #' @return a vector of characters.
  #'
  #' @import glue
  #'
  cell_ids <- c()
  counter <- list()
  for (cell_type in unique(labels)) {counter[[cell_type]] <- 1}
  for (n in 1:length(labels)) {
    cell_type <- labels[[n]]
    cell_ids[n] <- glue::glue("{cell_type}_{counter[[cell_type]]}")
    counter[[cell_type]] <- counter[[cell_type]] + 1}
  return(cell_ids)
}

get_data <- function(dataset) {
  #' Get a data.frame and the ground truth of a dataset loaded with TMExplorer.
  #'
  #' @param dataset a list loaded with TMExplorer::queryTME.
  #'
  #' @return a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @import SingleCellExperiment
  #' @import stats
  #' @import SummarizedExperiment
  #'
  throwaway_labels <- c("?", "NA", "", "Doublets")
  predictions <- SummarizedExperiment::colData(dataset)
  ground_truth <- stats::setNames(predictions$label, rownames(predictions))
  ground_truth <- ground_truth[!ground_truth %in% throwaway_labels]
  dataset <- SingleCellExperiment::counts(dataset)
  dataset <- dataset[, names(ground_truth)]
  data <- list(dataset=dataset, ground_truth=ground_truth)
  return(data)
}

sceve_LoadData <- function(dataset_label) {
  #' Load a scRNA-seq dataset of raw count expression.
  #'
  #' @param dataset_label one of `Li_HumCRC_b`, `Darmanis_HumGBM`, `JerbyArnon_HumMLM`,
  #' `Gillen_HumEPN`, `VanGalen_HumAML`, `Lambrechts_HumNSCLC` or `Peng_HumPDAC`.
  #' (cf. `sceve::get_datasets()`)
  #'
  #' @return a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #'
  #' @import TMExplorer
  #'
  #' @export
  #'
  datasets <- sceve_GetDatasets()
  accession_number <- datasets[dataset_label, "accession_number"]
  data <- TMExplorer::queryTME(accession_number)[[1]]

  if (accession_number == "GSE84465") {data <- head(data, -6)}
  # the 6 last rows are metadata in the Darmanis_HumGBM dataset.
  data <- get_data(data)

  cell_ids <- get_cell_ids(data$ground_truth)
  names(data$ground_truth) <- cell_ids
  data$ground_truth <- factor(data$ground_truth)

  colnames(data$dataset) <- cell_ids
  rownames(data$dataset) <- gsub("_", "+", rownames(data$dataset))
  return(data)
}
