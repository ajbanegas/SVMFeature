% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SVMFeature.R
\name{SVMFeature}
\alias{SVMFeature}
\title{SVMFeature Class (S3 Implementation)}
\usage{
SVMFeature(
  data,
  inputs,
  output,
  costs,
  pop_size,
  num_fea,
  n_iter = 10,
  max_time = 300,
  mode = "iters",
  objective = "distance-epsilon"
)
}
\arguments{
\item{data}{DataFrame. The input data.}

\item{inputs}{Character vector. The names of the input features.}

\item{output}{Character. The name of the output variable.}

\item{costs}{Numeric vector. The costs associated with the features.}

\item{pop_size}{Numeric. The population size.}

\item{num_fea}{Numeric. The number of features.}

\item{n_iter}{Numeric. The number of iterations (default: 10).}

\item{max_time}{Numeric. The maximum execution time (default: 300 seconds).}

\item{mode}{Character. "iters" for iteration-based execution, "time" for time-based execution.}

\item{objective}{Character. "distance-epsilon", "confusion-matrix", "distance-epsilon-costs" or "confusion-matrix-costs"}
}
\value{
An S3 object of class "SVMFeature".
}
\description{
Implements a feature selection algorithm using Support Vector Machines (SVM).
}
