#' Generate Metropolis-Hastings Sample For A Component of \eqn{\lambda}
#'
#' Generates a sample for a component of \eqn{\lambda} through Metropolis-Hastings
#' algorithm.
#'
#' @param x value from previous iteration (set to initial value for the first
#'          iteration)
#' @param m mean parameter from the conditional posterior of \eqn{\lambda}
#' @param s standard deviation parameter from the conditional posterior of
#'          \eqn{\lambda}
#' @param alpha type of skewness (1 or -1)
#'
#' @return value for the MCMC sample for \eqn{\lambda}
#' @export
#' @importFrom stats pnorm rnorm runif

metropolis = function(x, m, s, alpha){
  y = rnorm(1, mean = m, sd = s)
  if(runif(1)>(pnorm(alpha*y)/pnorm(alpha*x))) y = x
  return(y)
}
