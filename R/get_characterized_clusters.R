"Functions called to identify characterized clusters.

	2025/03/05 @yanisaspic"

add_leftover_cluster <- function(population, robust_clusters, selected_data, params) {
  #' Add a leftover cluster to the list of robust clusters.
  #'
  #' It corresponds to a group of samples unassigned to any robust cluster.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param robust_clusters list where every element is a robust pool of samples.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features` and `robustness`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a list where every element is a pool of samples.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features`, and `robustness`.
  #'
  #' @import glue
  #' @import stats
  #'
  n_clusters <- length(robust_clusters)
  samples_of_population <- colnames(selected_data$dataset)
  samples_in_robust_clusters <- unlist(sapply(robust_clusters, "[[", "samples"))
  leftover_samples <- setdiff(samples_of_population, samples_in_robust_clusters)
  if (length(leftover_samples) > 0) {
    leftover_cluster <- list(base_clusters=c(), clustering_methods=c(), robustness=0,
                             samples=leftover_samples, label=glue::glue("{population}.L"))
    leftover_cluster[["features"]] <- params$characteristic_features_strategy(leftover_cluster, selected_data, params)
    robust_clusters[[n_clusters + 1]] <- leftover_cluster}
  robust_clusters <- stats::setNames(robust_clusters, sapply(X=robust_clusters, FUN="[[", "label"))
  return(robust_clusters)
}

get_characterized_clusters <- function(population, robust_clusters, selected_data, params, figures) {
  #' Get characterized meta-clusters. They correspond to the sub-populations predicted by the clustering recursion.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param robust_clusters a list where every element is a pool of samples.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #' @param selected_data named list, with two names: `dataset` and `SeuratObject`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #'
  #' @return a list where every element is a characterized pool of samples.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features` and `robustness`.
  #'
  #' @import glue
  #' @import grDevices
  #'
  #' @export
  #'
  add_features <- function(cluster) {
    cluster[["features"]] <- params$characteristic_features_strategy(cluster, selected_data, params)
    return(cluster)}
  robust_clusters <- lapply(X=robust_clusters, FUN=add_features)

  characterized_clusters <- params$characterized_clusters_strategy(robust_clusters, selected_data, params)
  if (length(characterized_clusters) == 0) {return(list())}

  clusters <- add_leftover_cluster(population, characterized_clusters, selected_data, params)
  if (length(clusters) == 2) {
    characterized_clusters <- params$characterized_clusters_strategy(clusters, selected_data, params)
    if (length(characterized_clusters) < 2) {return(list())}}
  # if there is a single robust cluster and the leftovers, both must be differentially characterized to be counted.

  if (figures) {
    plot <- draw_characteristic_features(clusters)
    grDevices::pdf(file=glue::glue("{params$figures_path}/{population}_characteristic_features.pdf"))
    print(plot)
    grDevices::dev.off()}
  return(clusters)
}

generate_color_scale <- function(labels){
  #' Generate a vector of colors equal to the number of identities in the sample.
  #'
  #' This function is directly copied from the repository of the SCpubr package.
  #' c.f. https://github.com/enblacar/SCpubr/blob/main/R/utils.R
  #'
  #' @param labels a vector of cluster labels.
  #'
  #' @return a named vector of colors.
  #'
  #' @import colorspace
  #' @import grDevices
  #'

  # this section is added to have a cohesive colormap __________________________
  get_tail <- function(label) {
    tmp <- strsplit(label, split=".", fixed=TRUE)[[1]]
    tail <- tmp[length(tmp)]
    return(tail)}
  is_leftover <- function(label) {substr(label, nchar(label), nchar(label))=="L"}
  robust_labels <- labels[!sapply(X=labels, FUN=is_leftover)]
  n <- get_tail(robust_labels[length(robust_labels)])
  #_____________________________________________________________________________

  colors <- colorspace::qualitative_hcl(as.numeric(n), palette = "Dark 3")
  colors <- grDevices::col2rgb(colors)
  colors <- grDevices::rgb2hsv(colors)
  colors["v", ] <- colors["v", ] - 0.1
  colors["s", ] <- colors["s", ] + 0.2
  colors["s", ][colors["s", ] > 1] <- 1
  colors <- grDevices::hsv(h = colors["h", ],
                           s = colors["s", ],
                           v = colors["v", ],
                           alpha = 1)

  # this section is added to have a cohesive colormap __________________________
  get_index <- function(label) {as.numeric(get_tail(label))}
  robust_indexes <- sapply(X=robust_labels, FUN=get_index)
  colors <- colors[robust_indexes]
  if (is_leftover(labels[length(labels)])) {colors <- c(colors, "#404040FF")}
  #_____________________________________________________________________________

  names(colors) <- labels
  return(colors)
}

draw_characteristic_features <- function(characterized_clusters) {
  #' Get an upset-plot representing the characteristic features predicted in each meta-cluster.
  #'
  #' @param characterized_clusters a list where every element is a characterized pool of samples.
  #' The elements are named lists, with six names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label`, `features` and `robustness`.
  #'
  #' @return a plot.
  #'
  #' @import ggVennDiagram
  #' @import SCpubr
  #'
  #' @export
  #'
  get_features.cluster <- function(cluster) {names(cluster$features)}
  features <- lapply(X=characterized_clusters, FUN=get_features.cluster)
  has_elements <- function(feats) {length(feats) > 0}
  features <- Filter(f=has_elements, x=features)
  labels <- names(features)
  colormap <- generate_color_scale(labels)

  # these functions are used to improve the aesthetics of the plots_____________
  add_grid_colors <- function(plot) {
    grid_colors <- rep(NA, length(plot[[1]]$data$name))
    grid_colors[1:length(labels)] <- labels
    plot[[1]]$layers[[1]]$aes_params$colour <- NULL
    plot[[1]]$layers[[2]]$aes_params$colour <- NULL
    plot[[1]] <- plot[[1]] + ggplot2::aes(colour=grid_colors) + ggplot2::scale_colour_manual(values=colormap)
    plot[[1]] <- plot[[1]] + ggplot2::theme(legend.position="none")
    return(plot)}

  add_intersection_colors <- function(plot) {
    intersection_colors <- rep(NA, length(plot[[2]]$data$name))
    intersection_colors[1:length(labels)] <- labels
    plot[[2]] <- plot[[2]] + ggplot2::aes(fill=intersection_colors) + ggplot2::scale_fill_manual(values=colormap)
    plot[[2]] <- plot[[2]] + ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0, .05)))
    plot[[2]] <- plot[[2]] + ggplot2::theme(axis.text.y=ggplot2::element_text(vjust=0.25), legend.position="none")
    return(plot)}

  add_size_colors <- function(plot) {
    plot[[3]] <- plot[[3]] + ggplot2::aes(fill=labels) + ggplot2::scale_fill_manual(values=colormap)
    plot[[3]] <- plot[[3]] + ggplot2::theme(legend.position="none")
    return(plot)}
  #_____________________________________________________________________________

  plot <- ggVennDiagram::ggVennDiagram(features, category.names=labels, force_upset=TRUE,
                                       order.set.by="name", order.intersect.by="none",
                                       relative_height=1.6, relative_width=0.4)
  plot <- add_grid_colors(plot)
  plot <- add_intersection_colors(plot)
  plot <- add_size_colors(plot)
  return(plot)
}
