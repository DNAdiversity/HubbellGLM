#include <Rcpp.h>
using namespace Rcpp;

// Functions to calculate the mean and the variance of the polyseries link function
// associated with HubbellGLM
double polyseries_mean_single(int size, double alpha, double sigma){
  double mu = 1.0;
  for(int j = 1; j < size; j ++){
    mu += alpha/(alpha + std::pow(j, 1 - sigma));
  }
  return(mu);
}

double polyseries_var_single(int size, double alpha, double sigma){
  double v = 0.0;
  double term = 0.0;
  for(int j = 1; j < size; j ++){
    term = alpha/(alpha + std::pow(j, 1 - sigma));
    v += term * (1 - term);
  }
  return(v);
}

// [[Rcpp::export]]
NumericVector polyseries_mean(NumericVector size, NumericVector alpha, double sigma) {
  if (size.size() != alpha.size()) {
    // Throw an error if lengths are different
    stop("The lengths of 'size' and 'alpha' must be the same.");
  }
  int n  = size.size();
  NumericVector res(n);
  for (int j = 0; j < n; j++) {
    res[j] = polyseries_mean_single(size[j], alpha[j], sigma);
  }
  return res;
}

// [[Rcpp::export]]
NumericVector polyseries_var(NumericVector size, NumericVector alpha, double sigma) {
  if (size.size() != alpha.size()) {
    // Throw an error if lengths are different
    stop("The lengths of 'size' and 'alpha' must be the same.");
  }
  int n  = size.size();
  NumericVector res(n);
  for (int j = 0; j < n; j++) {
    res[j] = polyseries_var_single(size[j], alpha[j], sigma);
  }
  return res;
}


