% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Solution.R
\name{calculate_epsilon_objective}
\alias{calculate_epsilon_objective}
\title{Method to Calculate the Epsilon Objective}
\usage{
calculate_epsilon_objective(solution, data, output)
}
\arguments{
\item{solution}{Solution class object containing data and features for classification}

\item{data}{A data.frame with the input features.}

\item{output}{A vector with class labels.}
}
\value{
Solution class object with the epsilon objective calculated
}
\description{
This method calculates the epsilon objective by evaluating misclassified points.
}
