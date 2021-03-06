% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiomechanicsFunctions.R
\name{impulse}
\alias{impulse}
\title{Calculate impulse from ground reaction force data}
\usage{
impulse(time, GRF)
}
\arguments{
\item{\code{time}}{a vector of numerical data on the time sequence}

\item{\code{GRF}}{an array of columns for the force data (assumed that force data are already synchronized to the time data)}
}
\description{
Estimates impulse from data of force over time by calculating the area under the curve.
}
\details{
Impulse is a measure of the force applied over a specific time period. The time and force data should already be ordered so that the first row is the beginning of the trial and the last row is the end of the trial.
}
\examples{
time <- seq(1:10)
set.seed(123)
GRF <- data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))

impulse(time, GRF)


}
