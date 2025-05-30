% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/breve_GetCharacteristicFeatures.R
\name{get_differential_expression}
\alias{get_differential_expression}
\title{Get the features differentially expressed in an in-group (1),
with regards to an out-group (0).}
\usage{
get_differential_expression(filtered_dataset, groups)
}
\arguments{
\item{filtered_dataset}{an -omics dataset, with variant features.
Its rows are features and its columns are samples.}

\item{groups}{an ordered vector of integers.}
}
\value{
a data.frame associating genes to their log2 fold change and FDR.
}
\description{
Get the features differentially expressed in an in-group (1),
with regards to an out-group (0).
}
