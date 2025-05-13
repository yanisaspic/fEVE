"Functions used to get meta-clusters, i.e. robust clusters and a leftover cluster.

	2025/03/05 @yanisaspic"

get_transaction_database <- function(base_clusters) {
  #' Get a transaction database from the base clusters predicted.
  #'
  #' @param base_clusters a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @return a transactions object, where every sampleis a transaction and every predicted
  #' cluster is an item.
  #'
  #' @import arules
  #' @import glue
  #' @import rlang
  #' @import utils
  #'
  tmp <- rlang::hash(base_clusters) # use hash to prevent issues with parallel calculations.
  path <- glue::glue("./{tmp}.tmp")
  utils::write.table(base_clusters, file=path, col.names=FALSE, row.names=FALSE)
  transaction_database <- arules::read.transactions(path)
  file.remove(path)
  return(transaction_database)
}

get_associations <- function(transaction_database) {
  #' Measure the strength of the association between pairs of predicted clusters,
  #' according to the confidence metric of the frequent itemset mining framework.
  #'
  #' The confidence(A->C) corresponds to the proportion of samples in the cluster A, that
  #' are also in the cluster C. Note that confidence(A->C) can differ from confidence(C->A).
  #' Associations involving a minority of samples (i.e. confidence < 0.5) are filtered out.
  #'
  #' @param transaction_database a transactions object, where every sample is a transaction and every predicted
  #' cluster is an item.
  #'
  #' @return a data.frame associating clusters `A` and `C` to their `confidence(A->C)`.
  #'
  #' @import arules
  #'
  data <- arules::apriori(transaction_database,
                          parameter = list(support=0.001, confidence=0.5, minlen=2, maxlen=2, target="rules"))
  associations <- arules::DATAFRAME(data)

  get_item <- function(x) {substr(x, 2, nchar(x)-1)}
  associations <- list(A=sapply(X=as.character(associations$LHS), FUN=get_item),
                                C=sapply(X=as.character(associations$RHS), FUN=get_item),
                                confidence=as.numeric(associations$confidence))
  associations <- do.call(cbind, associations)
  associations <- data.frame(associations, row.names = NULL)
  return(associations)
}

get_strong_similarities <- function(associations) {
  #' Get pairs of clusters that share a strong similarity.
  #'
  #' A strong similarity is observed when two clusters share the majority of their samples,
  #' i.e. confidence(A->C) > 0.5 & confidence(C->A) > 0.5.
  #' It corresponds to the minimum proportion of samples shared between two clusters.
  #'
  #' @param associations a data.frame associating two clusters `A` and `C`, to their `confidence(A->C)`.
  #'
  #' @return a data.frame associating two clusters `A` and `C`, to their `similarity`.
  #'
  associations <- associations[order(associations$confidence, decreasing = TRUE),]

  get_id_row <- function(row) {
    elements <- sort(row[c("A", "C")])
    id <- paste(elements, collapse=".")
    return(id)}
  ids <- apply(X=associations, MARGIN=1, FUN=get_id_row)

  strong_similarities <- associations[duplicated(ids),]
  colnames(strong_similarities)[3] <- "similarity"
  return(strong_similarities)
}

get_samples_of_subgraph <- function(subgraph, base_clusters) {
  #' Get the samples at the intersection of all the clusters in a subgraph.
  #'
  #' @param subgraph a subgraph corresponding to strongly similar clusters.
  #' @param base_clusters a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #'
  #' @return a vector of samples.
  #'
  clusters_of_subgraph <- names(igraph::V(subgraph))
  cluster_is_in_subgraph <- function(cluster) {cluster %in% clusters_of_subgraph}
  data <- apply(X=base_clusters, MARGIN=c(1,2), FUN=cluster_is_in_subgraph)
  samples <- rownames(data[rowSums(data)==length(clusters_of_subgraph), ])
  return(samples)
}

get_subgraphs <- function(population, strong_similarities, base_clusters, params) {
  #' Model a graph where every node is a base cluster predicted, and edges are strong similarities.
  #'
  #' The graph is disjoint, and subgraphs (i.e. connected components) are extracted from it.
  #' They will be used to identify robust clusters.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param strong_similarities a data.frame associating two clusters `A` and `C`, to their `similarity`.
  #' @param base_clusters a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #'
  #' @return a list where every element is a pool of samples grouped together by multiple clustering methods.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #'
  #' @import igraph
  #' @import glue
  #'
  data <- igraph::graph.data.frame(strong_similarities)
  subgraphs <- igraph::decompose.graph(data)
  n_methods <- ncol(base_clusters)
  theoretical_robustness <- n_methods * (n_methods - 1) / 2

  # these functions are used to format the subgraph data in a list______________
  get_clusters_of_subgraph <- function(subgraph) {names(igraph::V(subgraph))}
  get_robustness_of_subgraph <- function(subgraph) {
    sum(as.numeric(igraph::E(subgraph)$similarity)) / theoretical_robustness}
  get_method_of_renamed_cluster <- function(renamed_cluster) {
    strsplit(renamed_cluster, split="_")[[1]][1]}
  get_methods_of_subgraph <- function(subgraph) {
    unique(sapply(X=get_clusters_of_subgraph(subgraph), FUN=get_method_of_renamed_cluster))}

  get_subgraph <- function(subgraph) {
    list(base_clusters=get_clusters_of_subgraph(subgraph),
         robustness=get_robustness_of_subgraph(subgraph),
         clustering_methods=get_methods_of_subgraph(subgraph),
         samples=get_samples_of_subgraph(subgraph, base_clusters))}
  #_____________________________________________________________________________

  subgraphs <- lapply(X=subgraphs, FUN=get_subgraph)
  subgraph_has_intersection <- function(subgraph) {length(subgraph$samples) > 0}
  subgraphs <- Filter(f=subgraph_has_intersection, x=subgraphs)
  robustnesses <- sapply(subgraphs, "[[", "robustness")
  subgraphs <- subgraphs[order(-robustnesses)]
  for (i in 1:length(subgraphs)) {subgraphs[[i]][["label"]] <- glue::glue("{population}.{i}")}
  return(subgraphs)
}

filter_conflicting_samples <- function(robust_clusters) {
  #' Filter out conflicting samples belonging to multiple robust clusters.
  #'
  #' If a sample exists in multiple robust clusters, it is exclusively associated
  #' to its most robust cluster. Robust clusters without samples are also removed.
  #'
  #' @param robust_clusters a list where every element is a pool of samples.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #'
  #' @return a list where every element is a pool of samples.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #'
  samples <- c()
  for (i in 1:length(robust_clusters)) {
    robust_clusters[[i]]$samples <- setdiff(robust_clusters[[i]]$samples, samples)
    samples <- c(samples, robust_clusters[[i]]$samples)}
  cluster_has_samples <- function(cluster) {length(cluster$samples) > 0}
  robust_clusters <- Filter(f=cluster_has_samples, x=robust_clusters)
  return(robust_clusters)
}

get_robustness_threshold <- function(population, base_clusters, records) {
  #' Get the robustness threshold to identify robust clusters.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param base_clusters a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #'
  #' @return a numeric.
  #'
  n_methods <- ncol(base_clusters)
  clusters_for_majority <- floor(n_methods / 2 + 1)
  edges_for_majority <- clusters_for_majority * (clusters_for_majority - 1) / 2
  majority_robustness <- 0.5 * edges_for_majority / (n_methods * (n_methods - 1) / 2)

  lineage_robustness <- records$meta[population, "robustness"]
  parent <- records$meta[population, "parent"]
  while (!is.na(parent)) {
    lineage_robustness <- max(records$meta[parent, "robustness"], lineage_robustness)
    parent <- records$meta[parent, "parent"]}

  robustness_threshold <- max(majority_robustness, lineage_robustness)
  return(robustness_threshold)
}

get_robust_clusters <- function(population, base_clusters, selected_data, records, params, figures) {
  #' Extract robust clusters and a leftover cluster from a set of base clusters predicted
  #' with multiple clustering methods.
  #'
  #' An empty list is returned if no robust cluster was predicted.
  #'
  #' @param population a character. It corresponds to the population that fEVE will attempt to cluster.
  #' @param base_clusters a data.frame associating samples to their predicted clusters.
  #' Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
  #' @param records a named list, with four data.frames: `samples`, `features`, `meta` and `methods`.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param figures a boolean that indicates if figures should be drawn to explain the clustering recursion.
  #'
  #' @return a list where every element is a pool of samples.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #'
  #' @import glue
  #' @import grDevices
  #' @import stats
  #'
  #' @export
  #'
  transaction_database <- get_transaction_database(base_clusters)
  associations <- get_associations(transaction_database)
  strong_similarities <- get_strong_similarities(associations)
  subgraphs <- get_subgraphs(population, strong_similarities, base_clusters, params)

  robustness_threshold <- get_robustness_threshold(population, base_clusters, records)
  subgraph_is_robust_cluster <- function(subgraph) {subgraph$robustness > robustness_threshold}
  robust_clusters <- Filter(f=subgraph_is_robust_cluster, x=subgraphs)
  robust_clusters <- stats::setNames(robust_clusters, sapply(X=robust_clusters, FUN="[[", "label"))

  if (length(robust_clusters) > 0 & figures) {
    plot <- draw_robust_clusters(robust_clusters, selected_data)
    grDevices::pdf(file=glue::glue("{params$figures_path}/{population}_robust_clusters.pdf"))
    print(plot)
    grDevices::dev.off()}
  return(robust_clusters)
}

draw_robust_clusters <- function(robust_clusters, selected_data) {
  #' Get composite U-MAP plots representing the meta-clusters predicted with the fEVE framework.
  #'
  #' @param robust_clusters a list where every element is a pool of samples.
  #' The elements are named lists, with five names:
  #' `base_clusters`, `samples`, `clustering_methods`, `label` and `robustness`.
  #' @param selected_data a named list, with two names: `dataset` and `SeuratObject`.
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
  get_colormap_cluster <- function(cluster) {
    plot_label <- glue::glue("{cluster$label} ({round(cluster$robustness, 2)})")
    stats::setNames(rep(plot_label, length(cluster$samples)), cluster$samples)}
  plot_coloring <- lapply(X=robust_clusters, FUN=get_colormap_cluster)
  plot_coloring <- unlist(unname(plot_coloring))
  selected_data$SeuratObject[["robust_clusters"]] <- as.factor(plot_coloring)

  plot <- SCpubr::do_DimPlot(selected_data$SeuratObject,
                             split.by="robust_clusters",
                             legend.position="none",
                             na.value="grey25")
  for (i in 1:length(plot)) {
    plot[[i]][[1]] <- plot[[i]][[1]] +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, margin=ggplot2::margin(1, 0, 0, 0)),
                     panel.background=ggplot2::element_rect(fill="lightgrey"),
                     legend.position="none", axis.title=ggplot2::element_blank())}
  return(plot)
}
