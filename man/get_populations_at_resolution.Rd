% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_populations_at_resolution}
\alias{get_populations_at_resolution}
\title{Get all the populations at a specific resolution.}
\usage{
get_populations_at_resolution(resolution, records_samples)
}
\arguments{
\item{resolution}{an integer.}

\item{records_samples}{a data.frame associating samples to their predicted populations.
Its rows are samples and its columns are population. The cell values range from 0 to 1.}
}
\value{
a vector of population labels.
}
\description{
The root population is resolution 1, and its children populations are resolution 2, etc.
}
