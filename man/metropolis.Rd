% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metropolis.R
\name{metropolis}
\alias{metropolis}
\title{Generate Metropolis-Hastings Sample For A Component of \eqn{\lambda}}
\usage{
metropolis(x, m, s, alpha)
}
\arguments{
\item{x}{value from previous iteration (set to initial value for the first
iteration)}

\item{m}{mean parameter from the conditional posterior of \eqn{\lambda}}

\item{s}{standard deviation parameter from the conditional posterior of
\eqn{\lambda}}

\item{alpha}{type of skewness (1 or -1)}
}
\value{
value for the MCMC sample for \eqn{\lambda}
}
\description{
Generates a sample for a component of \eqn{\lambda} through Metropolis-Hastings
algorithm.
}
