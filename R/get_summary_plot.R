"Functions called to draw the summary figure of a clustering analysis.

	2025/02/24 @yanisaspic"

get_tree_data <- function(records_meta) {
  #' Get a data structure representing the hierarchy of clusters predicted.
  #' The robustness of each cluster is also included.
  #' The resulting data structure is suitable for `ggtree::ggtree()`.
  #' 
  #' @param records_meta associating predicted populations to generic information, including:
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

get_tree_plot <- function(tree_data) {
  #' Get a plot representing the hierarchy of clusters predicted.
  #' 
  #' @param tree_data a data structure suitable for `ggtree::ggtree()`.
  #'
  #' @return a plot.
  #' 
  #' @import ggplot2
  #' @import ggtree
  #' 
  is_robust <- function(robustness) {ifelse(robustness > 0, "yes", "no")}
  tree_data[, "is_robust"] <- sapply(X=tree_data$robustness, FUN=is_robust)
  
  plot <- ggtree::ggtree(tree_data) +
    ggplot2::geom_label(ggplot2::aes(label=label, fill=is_robust, color=is_robust),
                        hjust=ifelse(tree_data$isTip, "right", "middle")) +
    ggplot2::scale_color_manual(values=c("white", "black")) +
    ggplot2::scale_fill_manual(values=list("yes"="grey90", "no"="grey10")) +
    ggplot2::geom_label(data=tree_data[-1,], ggplot2::aes(x=branch, label=round(robustness, 2)),
                        label.size=NA, hjust="right") +
    ggplot2::theme(legend.position="none")
  
  return(plot)
}

get_cluster_sizes <- function(records_samples) {
  #' Get a data.frame associating each cluster to its size.
  #' 
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and and its columns are population. The cell values range from 0 to 1.
  #' 
  #' @return a data.frame with two columns: `cluster` and `n`.
  #'   
  get_cluster_size <- function(cluster) {
    tmp <- get_samples_of_population(cluster, records_samples)
    cluster_size <- data.frame(cluster=cluster, n=length(tmp))
    return(cluster_size)}
  
  cluster_sizes <- lapply(X=colnames(records_samples), FUN=get_cluster_size)
  cluster_sizes <- do.call(rbind, cluster_sizes)
  return(cluster_sizes)
}

get_barplot_sizes <- function(tree_data, records_samples) {
  #' Get a barplot associating each leaf cluster to its size.
  #' 
  #' @param tree_data a data structure suitable for `ggtree::ggtree()`.
  #' @param records_samples a data.frame associating samples to their predicted populations.
  #' Its rows are samples and and its columns are population. The cell values range from 0 to 1.
  #' 
  #' @return a plot.
  #' 
  #' @import ggplot2
  #' 
  leaves <- tree_data[tree_data$isTip, ]
  alignment <- leaves$label[order(leaves$y)]
  cluster_sizes <- get_cluster_sizes(records_samples)
  cluster_sizes <- cluster_sizes[cluster_sizes$cluster %in% leaves$label, ]
  
  plot <- ggplot2::ggplot(data=cluster_sizes) +
    ggplot2::geom_bar(ggplot2::aes(x=factor(cluster, levels=alignment), y=n),
                      stat="identity", fill="black") +
    ggplot2::geom_text(aes(x=factor(cluster, levels=alignment), y=n, label=n), color="white",
                       position = ggplot2::position_stack(0.5), fontface="bold") +
    ggplot2::scale_y_log10(expand=ggplot2::expansion(mult=0), guide="axis_logticks") +
    ggplot2::coord_flip()

  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::ylab("# samples")
  return(plot)
}

get_cluster_features <- function(records_features) {
  #' Get a data.frame associating each cluster to its features.
  #' 
  #' @param records_features a data.frame associating predicted populations to their characteristic features.
  #' Its rows are genes, its columns are predicted populations,
  #' and the strength of the characterization (e.g. log2-transformed pvalues) are reported in the table.
  #' 
  #' @return a data.frame with three columns: `cluster`, `n` and `names`.
  #'   
  get_cluster_features <- function(cluster) {
    tmp <- records_features[records_features[, cluster] > 0,]
    if (nrow(tmp) == 0) {
        cluster_features <- data.frame(cluster=cluster, n=0, label=NA)
        return(cluster_features)}
    tmp <- tmp[order(tmp[, cluster], decreasing=TRUE), ]
    cluster_features <- data.frame(cluster=cluster, n=nrow(tmp), label=rownames(tmp)[1])
    return(cluster_features)}
  
  cluster_features <- lapply(X=colnames(records_features), FUN=get_cluster_features)
  cluster_features <- do.call(rbind, cluster_features)
  return(cluster_features)
}

get_barplot_features <- function(tree_data, records_features) {
  #' Get a barplot associating each leaf cluster to its features.
  #' 
  #' @param tree_data a data structure suitable for `ggtree::ggtree()`.
  #' @param records_features a data.frame associating predicted populations to their characteristic features.
  #' Its rows are genes, its columns are predicted populations,
  #' and the strength of the characterization (e.g. log2-transformed pvalues) are reported in the table.
  #' 
  #' @return a plot.
  #' 
  #' @import ggplot2
  #' 
  leaves <- tree_data[tree_data$isTip, ]
  alignment <- leaves$label[order(leaves$y)]
  cluster_features <- get_cluster_features(records_features)
  cluster_features <- cluster_features[cluster_features$cluster %in% leaves$label, ]
  
  plot <- ggplot2::ggplot(data=cluster_features) +
    ggplot2::geom_bar(ggplot2::aes(x=factor(cluster, levels=alignment), y=n),
                      stat="identity", fill="black") +
    ggplot2::geom_text(aes(x=factor(cluster, levels=alignment), y=n, label=label), color="white",
                       position = ggplot2::position_stack(0.5), fontface="bold") +
    ggplot2::scale_y_log10(expand=ggplot2::expansion(mult=0), guide="axis_logticks") +
    ggplot2::coord_flip()

  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::ylab("# features")
  return(plot)
}

get_summary_plot <- function(records, widths=c(8, 4, 4)) {
  #' Get a composite plot summarizing a fEVE clustering analysis.
  #' The composite reports the relationships, as well as the sizes and the markers of the leaf clusters.
  #' 
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
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
  tree_plot <- get_tree_plot(tree_data)
  tree_plot <- tree_plot +
    ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0, 0, 0), "null")) +
    ggplot2::guides(fill="none", color="none")
  
  remove_y_axis <- function(barplot) {
    barplot <- barplot +
      ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank())
    return(barplot)}
  
  barplot_sizes <- get_barplot_sizes(tree_data, records$samples)
  barplot_features <- get_barplot_features(tree_data, records$features)
  
  composite_plot <- tree_plot +
    remove_y_axis(barplot_sizes) +
    remove_y_axis(barplot_features) +
    patchwork::plot_layout(widths=widths, guides="collect") &
    ggplot2::theme(legend.position="bottom")
  
  composite_plot <- composite_plot +
      patchwork::plot_annotation(tag_levels="a")
  return(composite_plot)
}