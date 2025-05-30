% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_characterized_clusters.R
\name{add_leftover_cluster}
\alias{add_leftover_cluster}
\title{Add a leftover cluster to the list of robust clusters.}
\usage{
add_leftover_cluster(population, robust_clusters, selected_data, params)
}
\arguments{
\item{population}{a character. It corresponds to the population that fEVE will attempt to cluster.}

\item{robust_clusters}{list where every element is a robust pool of samples.
The elements are named lists, with six names:
\code{base_clusters}, \code{samples}, \code{clustering_methods}, \code{label}, \code{features} and \code{robustness}.}

\item{selected_data}{a named list, with two names: \code{dataset} and \code{SeuratObject}.}

\item{params}{a list of parameters (cf. \code{feve::get_parameters()}).}
}
\value{
a list where every element is a pool of samples.
The elements are named lists, with six names:
\code{base_clusters}, \code{samples}, \code{clustering_methods}, \code{label}, \code{features}, and \code{robustness}.
}
\description{
It corresponds to a group of samples unassigned to any robust cluster.
}
