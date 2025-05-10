"Function called to plot the contributions of each clustering method in identifying robust clusters.

	2025/05/09 @yanisaspic"

get_barplot_robustness <- function(records_meta) {
  #' Get a barplot associating each robust cluster to its robustness value.
  #'
  #' @param records_meta a data.frame associating predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #'
  #' @return a plot.
  #'
  #' @import ggplot2
  #'
  data <- records_meta[records_meta$robustness > 0,]
  data[, "cluster"] <- rownames(data)

  plot <- ggplot2::ggplot(data) +
    ggplot2::geom_bar(ggplot2::aes(x=cluster, y=robustness, fill=robustness), stat="identity") +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=0)) +
    ggplot2::scale_fill_gradient(low="#EBEBEB", high="#0072B2", limits=c(0,1))

  plot <- plot + ggplot2::theme_classic() +
    ggplot2::theme(legend.position="none",
                   axis.text.y=ggplot2::element_text(vjust=0))
  return(plot)
}

get_heatmap_contributions <- function(records_methods) {
  #' Get a heatmap plot representing which method contributed to the prediction of each cluster.
  #'
  #' @param records_methods a data.frame associating methods (rows) to the robust clusters
  #' they predicted (columns). Cell values are binary, and 1 indicates that a method contributed
  #' in the detection of a robust cluster.
  #'
  #' @return a plot
  #'
  #' @import ggplot2
  #' @import tibble
  #' @import tidyr
  #'
  data <- tibble::rownames_to_column(records_methods, var="method")
  data <- tidyr::pivot_longer(data, cols=-method, names_to="cluster", values_to="value")

  plot <- ggplot2::ggplot(data) +
    ggplot2::geom_tile(ggplot2::aes(x=cluster, y=method, fill=value), color="white", linewidth=1) +
    ggplot2::scale_fill_gradient(low="#EBEBEB", high="#0072B2", limits=c(0,1)) +
    ggplot2::scale_x_discrete(expand=ggplot2::expansion(mult=0.11)) +
    ggplot2::scale_y_discrete(expand=ggplot2::expansion(mult=0.18))

  plot <- plot + ggplot2::theme_classic() +
    ggplot2::theme(legend.position="none",
                   panel.border = ggplot2::element_rect(colour="black", fill=NA, linewidth=1),
                   axis.text.x = ggplot2::element_text(angle=45, vjust=1, hjust=1),
                   axis.title.x = ggplot2::element_blank())
  return(plot)
}

get_plot_contributions <- function(records, heights=c(1,1)) {
  #' Get a composite plot representing the contribution of each clustering method
  #' on the fEVE prediction.
  #'
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param heights a vector of numeric.
  #'
  #' @return a composite plot
  #'
  #' @import ggplot2
  #' @import patchwork
  #'
  #' @export
  #'
  subplots <- list(
    "contributions"=get_heatmap_contributions(records$methods),
    "robustness"=get_barplot_robustness(records$meta))

  remove_margins <- function(plot) {
    plot <- plot + ggplot2::theme(plot.margin=ggplot2::unit(c(0,0,0,0), "cm"))}
  for (p in names(subplots)) {subplots[[p]] <- remove_margins(subplots[[p]])}

  remove_x_axis <- function(plot) {
    plot <- plot +
      ggplot2::theme(axis.line.x=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank())}

  composite_plot <- remove_x_axis(subplots$robustness) + subplots$contributions +
    patchwork::plot_layout(ncol=1, nrow=2, guides="collect", heights=heights)
  return(composite_plot)
}
