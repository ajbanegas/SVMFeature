% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Population.R
\name{generate_initial_population}
\alias{generate_initial_population}
\title{Generate initial population}
\usage{
generate_initial_population(population, data, inputs, output, costs)
}
\arguments{
\item{population}{Population class object}

\item{data}{A data.frame with the input features}

\item{inputs}{A vector of input features}

\item{output}{A vector with class labels}

\item{costs}{A vector with feature costs}
}
\value{
Population class object updated with the generated initial population
}
\description{
This function generates the initial population.
}
