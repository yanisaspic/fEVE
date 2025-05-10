"Function used to generate a report file with summary figures

	2025/05/09 @yanisaspic"

get_plot_wrapper <- function(plot_label, data, records) {
  #' Wrapper function to generate and save a plot.
  #'
  #' @param plot_label a character.
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a dataset, without selected features. Its rows are features and its columns are samples.
  #' `ground_truth` is a named factor associating samples to their cluster annotations.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @return a named list with two names: `plot` and `dimensions`.
  #'
  if (plot_label == "summary") {plot <- get_plot_summary(records)}
  if (plot_label == "compositions") {plot <- get_plot_compositions(records, data$ground_truth)}
  if (plot_label == "contributions") {plot <- get_plot_contributions(records)}
  if (plot_label == "resolutions") {plot <- get_plot_resolutions(data, records$samples)}

  if (plot_label %in% c("summary", "compositions")) {dimensions <- list(width=16, height=9)}
  if (plot_label %in% c("contributions", "resolutions")) {dimensions <- list(width=3.5, height=4)}
  out <- list(plot=plot, dim=dimensions)
  return(out)
}

get_placeholder_ground_truth <- function(records, random_state=1) {
  #' Get a placeholder ground truth corresponding to shuffled fEVE labels.
  #'
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param random_state a numeric.
  #'
  #' @return a named factor associating samples to their ground truth label.
  #'
  #' @import stats
  #'
  set.seed(random_state)
  preds <- get_leaf_clusters(records$samples)
  labels <- sample(as.vector(preds), length(preds))
  ground_truth <- stats::setNames(labels, names(preds))
  return(ground_truth)
}

get_analysis_report <- function(dataset_init, records, ground_truth=NULL, path="./report.pdf") {
  #' Save summary figures of an analysis in a .pdf file.
  #'
  #' @param dataset_init a dataset, without selected features.
  #' Its rows are features and its columns are samples.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param ground_truth a named factor associating samples to their ground truth label.
  #' @param path a character.
  #'
  #' @import glue
  #' @import grid
  #' @import grDevices
  #' @import png
  #' @import stats
  #'
  #' @export
  #'
  if (is.null(ground_truth)) {ground_truth <- get_placeholder_ground_truth(records)}
  data <- list(dataset=dataset_init, ground_truth=ground_truth)

  draw_plot_png <- function(plot_label) {
    tmp <- get_plot_wrapper(plot_label, data, records)
    plot_path <- glue::glue("./{plot_label}.png")
    ggplot2::ggsave(plot_path, tmp$plot, width=tmp$dim$width, height=tmp$dim$height, dpi=300)
    return(plot_path)}
  plot_labels <- c("summary", "compositions", "contributions", "resolutions")
  drawings_paths <- sapply(X=plot_labels, FUN=draw_plot_png)

  draw_plot_pdf <- function(drawing_path) {
    img <- png::readPNG(drawing_path)
    grid::grid.newpage()
    grid::grid.raster(img)}
  grDevices::pdf(file=path)
  sapply(X=drawings_paths, FUN=draw_plot_pdf)
  grDevices::dev.off()
  unlink(drawings_paths)
}
