"Function called to plot the cluster quality across multiple resolutions, with multiple metrics.

	2025/05/09 @yanisaspic"

get_resolution_data <- function(data, records_samples, method_label) {
  #' Get the clustering metrics at each resolution to be plotted.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a dataset, without selected features. Its rows are features and its columns are samples.
  #' `ground_truth` is a named factor associating samples to their cluster annotations.
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are population. The cell values range from 0 to 1.
  #' @param method_label a character.
  #'
  #' @return a data.frame with five columns: `method`, `n_samples`, `resolution`, `metric` and `value`.
  #'
  #' @import tidyr
  #'
  f <- function(resolution) {
    preds <- get_leaf_clusters_at_resolution(resolution, records_samples)
    benchmark <- get_clustering_metrics_trifecta(data, preds, method_label)
    benchmark[, "resolution"] <- resolution
    return(benchmark)}

  maximum_resolution <- get_maximum_resolution(records_samples)
  results <- lapply(X=1:maximum_resolution, FUN=f)
  results <- do.call(rbind, results)
  results <- tidyr::pivot_longer(results, cols=c("ARI", "NMI", "nPurity", "SI"), names_to="metric")
  results$value <- as.numeric(results$value)
  return(results)
}

get_plot_resolutions <- function(data, records_samples, method_label="fEVE") {
  #' Get lineplots summarizing the cluster quality of a fEVE prediction across multiple resolutions.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a dataset, without selected features. Its rows are features and its columns are samples.
  #' `ground_truth` is a named factor associating samples to their cluster annotations.
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are population. The cell values range from 0 to 1.
  #' @param method_label a character.
  #'
  #' @return a faceted plot.
  #'
  #' @import ggplot2
  #'
  #' @export
  #'
  resolution_data <- get_resolution_data(data, records_samples, method_label)

  plot <- ggplot2::ggplot(resolution_data, ggplot2::aes(x=resolution, y=value)) +
    ggplot2::geom_line(ggplot2::aes(group=method, color=method)) +
    ggplot2::geom_point(ggplot2::aes(group=method, color=method)) +
    ggplot2::facet_wrap(metric ~ ., scales="free_y") +
    ggplot2::scale_color_brewer(palette="Dark2") +
    ggplot2::scale_x_continuous(breaks=seq(1, max(resolution_data$resolution), by=1))

  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid.major=ggplot2::element_line(linewidth=0.5),
                   panel.grid.minor=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(), legend.text=ggplot2::element_text(size=10),
                   legend.title=ggplot2::element_blank(), legend.position="top",
                   strip.background=ggplot2::element_blank(),
                   panel.border=ggplot2::element_rect(colour="black", fill=NA))
  return(plot)
}
