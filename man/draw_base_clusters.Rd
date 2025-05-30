% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_base_clusters.R
\name{draw_base_clusters}
\alias{draw_base_clusters}
\title{Get composite U-MAP plots representing the base clusters predicted by each clustering method.}
\usage{
draw_base_clusters(selected_SeuratObject, base_clusters)
}
\arguments{
\item{selected_SeuratObject}{a SeuratObject, on which the function ScaleData()
of Seurat has been applied already.}

\item{base_clusters}{a data.frame associating samples to their predicted clusters.
Its rows are samples, its columns are clustering methods, and predicted populations are reported in the table.}
}
\value{
a plot.
}
\description{
Get composite U-MAP plots representing the base clusters predicted by each clustering method.
}
