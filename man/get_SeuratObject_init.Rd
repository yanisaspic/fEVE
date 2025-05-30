% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_SeuratObject_init}
\alias{get_SeuratObject_init}
\title{Get a SeuratObject from a dataset, without selected features.}
\usage{
get_SeuratObject_init(dataset_init)
}
\arguments{
\item{dataset_init}{a dataset, without selected features.
Its rows are features and its columns are samples.}
}
\value{
a SeuratObject, on which the function RunUMAP() of Seurat has been applied already.
}
\description{
This function is used once prior to a fEVE clustering analysis in order to draw the
extracted data at each clustering recursion (cf. \code{feve::draw_extracted_data()}).
}
