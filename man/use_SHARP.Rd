% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sceve_GetBaseClusters.R
\name{use_SHARP}
\alias{use_SHARP}
\title{Predict clusters with the SHARP method.}
\usage{
use_SHARP(selected_dataset, params)
}
\arguments{
\item{selected_dataset}{a scRNA-seq dataset of raw count expression, with selected genes.}

\item{params}{a list of parameters (cf. \code{feve::get_parameters()}).}
}
\value{
a named factor that associates each cell to their cluster prediction.
}
\description{
Predict clusters with the SHARP method.
}
