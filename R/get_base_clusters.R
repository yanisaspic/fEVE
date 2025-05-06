"Functions used to get base clusters.

	2025/03/05 @yanisaspic"

draw_base_clusters <- function(selected_SeuratObject, base_clusters) {
  #' Get composite U-MAP plots representing the base clusters predicted by each clustering method.
  #'
  #' @param selected_SeuratObject a SeuratObject, on which the function ScaleData()
  #' of Seurat has been applied already.
  #' @param base_clusters a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @return a plot.
  #'
  #' @import ggplot2
  #' @import glue
  #' @import gridExtra
  #' @import SCpubr
  #' @import stats
  #'
  #' @export
  #'
  get_number <- function(renamed_cluster) {strsplit(renamed_cluster, split="_")[[1]][2]}
  base_clusters <- apply(base_clusters, c(1,2), get_number)

  get_plot_method <- function(method) {
    samples <- rownames(base_clusters)
    numbers <- factor(base_clusters[, method])
    selected_SeuratObject@active.ident <- stats::setNames(numbers, samples)
    n_clusters <- length(levels(numbers))
    plot <- SCpubr::do_DimPlot(selected_SeuratObject) +
      ggplot2::ggtitle(glue::glue("{method}: {n_clusters} cluster(s)")) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, margin=ggplot2::margin(1, 0, 0, 0)),
                     panel.background=ggplot2::element_rect(fill="lightgrey"), legend.position="none",
                     axis.title=ggplot2::element_blank())
    return(plot)
  }

  clustering_methods <- colnames(base_clusters)
  plots <- lapply(X=clustering_methods, FUN=get_plot_method)
  names(plots) <- clustering_methods
  composite_plot <- do.call(gridExtra::grid.arrange, plots)
  return(composite_plot)
}

get_base_clusters <- function(population, selected_data, params, figures) {
  #' Get base clusters predicted with multiple clustering methods.
  #'
  #' @param population a character. It corresponds to the population that feature will attempt to cluster.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #'
  #' @return a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @import glue
  #' @import grDevices
  #'
  #' @export
  #'
  base_clusters <- params$base_clusters_strategy(selected_data, params)
  if (figures) {
    grDevices::pdf(file=glue::glue("{params$figures_path}/{population}_base_clusters.pdf"))
    plot <- draw_base_clusters(selected_data$SeuratObject, base_clusters)
    print(plot)
    grDevices::dev.off()
  }
  return(base_clusters)
}
