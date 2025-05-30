% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_pending_population}
\alias{get_pending_population}
\title{Get a population for which no fEVE clustering recursion has been attempted.}
\usage{
get_pending_population(records_meta)
}
\arguments{
\item{records_meta}{a data.frame associating predicted populations to generic information, including:
their \code{size}, their \code{robustness}, their \code{parent} and their \code{clustering_status}.}
}
\value{
a character.
}
\description{
Get a population for which no fEVE clustering recursion has been attempted.
}
