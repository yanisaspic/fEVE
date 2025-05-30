% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_characterized_clusters.R
\name{get_characterized_clusters}
\alias{get_characterized_clusters}
\title{Get characterized meta-clusters. They correspond to the sub-populations predicted by the clustering recursion.}
\usage{
get_characterized_clusters(
  population,
  robust_clusters,
  selected_data,
  params,
  figures
)
}
\arguments{
\item{population}{a character. It corresponds to the population that fEVE will attempt to cluster.}

\item{robust_clusters}{a list where every element is a pool of samples.
The elements are named lists, with five names:
\code{base_clusters}, \code{samples}, \code{clustering_methods}, \code{label} and \code{robustness}.}

\item{selected_data}{named list, with two names: \code{dataset} and \code{SeuratObject}.}

\item{params}{a list of parameters (cf. \code{feve::get_parameters()}).}

\item{figures}{a boolean that indicates if figures should be drawn to explain the clustering recursion.}
}
\value{
a list where every element is a characterized pool of samples.
The elements are named lists, with six names:
\code{base_clusters}, \code{samples}, \code{clustering_methods}, \code{label}, \code{features} and \code{robustness}.
}
\description{
Get characterized meta-clusters. They correspond to the sub-populations predicted by the clustering recursion.
}
