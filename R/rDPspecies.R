# Function to simulate from the numer of distinct species from a Dirichlet process
#' @export
#' @importFrom Rcpp sourceCpp
#' @useDynLib HubbellGLM
rDPspecies <- Vectorize(FUN = DPspecies_cpp,
                        vectorize.args = c("alpha", "size"))
