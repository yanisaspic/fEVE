"Function called to draw the summary figure of a clustering analysis,
with regards to the ground truth of the dataset..

	2025/05/09 @yanisaspic"

get_tree_data <- function(records_meta) {
  #' Get a data structure representing the hierarchy of clusters predicted.
  #' The robustness of each cluster is also included.
  #' The resulting data structure is suitable for `ggtree::ggtree()`.
  #'
  #' @param records_meta a data.frame associating predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #'
  #' @return a data structure suitable for `ggtree::ggtree()`.
  #'
  #' @import dplyr
  #' @import ggplot2
  #' @import ggtree
  #' @import tidytree
  #'
  relationships <- data.frame(parent=records_meta[-1,]$parent, node=rownames(records_meta[-1,]))
  tmp <- tidytree::as.phylo(relationships)
  tree_data <- dplyr::as_tibble(tmp)
  coordinates <- ggplot2::fortify(tmp)

  # align the node clusters to their specific rows______________________________
  alignment <- rownames(records_meta)[order(records_meta$robustness, -records_meta$size)]
  records_meta <- records_meta[alignment, ]
  tree_data <- tree_data[match(alignment, tree_data$label), ]
  coordinates <- coordinates[match(alignment, coordinates$label), ]
  # ____________________________________________________________________________

  tree_data[, "robustness"] <- records_meta$robustness
  for (col in c("isTip", "x", "y", "branch", "angle")) {tree_data[, col] <- coordinates[, col]}
  return(tree_data)
}

get_cluster_compositions <- function(ground_truth, records_samples) {
  #' Get a data.frame associating each cluster to its class composition,
  #' according to some ground truth.
  #'
  #' @param ground_truth a named factor associating samples to their ground truth label.
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are populations. The cell values range from 0 to 1.
  #'
  #' @return a data.frame with three columns: `cluster`, `sample_class` and `n`.
  #'
  get_ground_truth <- function(sample_id) {ground_truth[sample_id]}

  get_cluster_comp <- function(cluster) {
    samples_of_cluster <- get_samples_of_population(cluster, records_samples)
    labels <- sapply(X=samples_of_cluster, FUN=get_ground_truth)
    tmp <- table(labels)
    cluster_comp <- data.frame(cluster=cluster, sample_type=names(tmp), n=as.numeric(tmp))
    return(cluster_comp)}

  cluster_compositions <- lapply(X=colnames(records_samples), FUN=get_cluster_comp)
  cluster_compositions <- do.call(rbind, cluster_compositions)
  return(cluster_compositions)
}

get_barplot_compositions <- function(tree_data, ground_truth, records_samples) {
  #' Get a barplot associating each leaf cluster to its sample composition.
  #'
  #' @param tree_data a data structure suitable for `ggtree::ggtree()`.
  #' @param ground_truth a named factor associating samples to their ground truth label.
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and its columns are populations. The cell values range from 0 to 1.
  #'
  #' @return a plot.
  #'
  #' @import ggplot2
  #' @import pals
  #' @import scales
  #'
  leaves <- tree_data[tree_data$isTip, ]
  alignment <- leaves$label[order(leaves$y)]
  cluster_compositions <- get_cluster_compositions(ground_truth, records_samples)
  cluster_compositions <- cluster_compositions[cluster_compositions$cluster %in% leaves$label, ]
  n_sample_types <- length(unique(cluster_compositions$sample_type))

  plot <- ggplot2::ggplot(data=cluster_compositions) +
    ggplot2::geom_bar(ggplot2::aes(x=factor(cluster, levels=alignment), y=n, fill=sample_type),
                      position="fill", stat="identity", color="black") +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=0), labels=scales::percent) +
    ggplot2::scale_fill_manual(values=pals::kelly(n_sample_types)) +
    ggplot2::coord_flip()

  n_rows_legend <- ceiling(n_sample_types / 7)
  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::ylab("sample types") +
    ggplot2::guides(fill=ggplot2::guide_legend(nrow=n_rows_legend))
  return(plot)
}

get_plot_compositions <- function(records, ground_truth, widths=c(8, 4, 4)) {
  #' Get a composite plot summarizing a scEVE clustering analysis on a dataset with a ground truth.
  #' The composite reports the relationships, the cell composition and the size of every leaf cluster
  #' predicted by scEVE.
  #'
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param ground_truth a named factor associating samples to their class.
  #' @param widths a vector of numeric.
  #'
  #' @return a composite plot.
  #'
  #' @import ggplot2
  #' @import patchwork
  #'
  #' @export
  #'
  tree_data <- get_tree_data(records$meta)
  plot_ground_truth <- get_plot_summary(records, widths)
  barplot_compositions <- get_barplot_compositions(tree_data, ground_truth, records$samples)

  remove_y_axis <- function(barplot) {
    barplot <- barplot +
      ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank())
    return(barplot)}
  barplot_compositions <- remove_y_axis(barplot_compositions)

  plot_ground_truth[[3]] <- plot_ground_truth[[2]]
  plot_ground_truth[[2]] <- barplot_compositions
  return(plot_ground_truth)
}
