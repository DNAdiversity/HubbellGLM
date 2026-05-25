#' Simulate the number of distinct species from a Dirichlet process
#'
#' Draws the number of distinct species \eqn{Y^{(n)}} from the Dirichlet
#' process prior with concentration parameter \eqn{\alpha} and sample size
#' \eqn{n}.
#'
#' @param alpha Concentration parameter (positive scalar or vector).
#' @param size Sample size \eqn{n} (positive integer, scalar or vector).
#'
#' @return An integer vector of simulated species counts.
#'
#' @export
#' @importFrom Rcpp sourceCpp
#' @useDynLib HubbellGLM
rDPspecies <- Vectorize(FUN = DPspecies_cpp,
                        vectorize.args = c("alpha", "size"))
