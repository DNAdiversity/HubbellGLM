#include <Rcpp.h>
using namespace Rcpp;

// Functions to calculate the mean and the variance of the polyseries link function
// associated with HubbellGLM
// [[Rcpp::export]]
double polyseries_mean_single(int size, double alpha, double sigma){
  double mu = 1.0;
  for(int j = 1; j < size; j ++){
    mu += alpha/(alpha + std::pow(j, 1 - sigma));
  }
  return(mu);
}

// [[Rcpp::export]]
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
NumericVector polyseries_mean_grid(NumericVector points_grid, int size, double alpha, double sigma) {
  int tot = points_grid.size();
  int id = 0;
  NumericVector res(tot);
  double mu = 0;
  for(int j = 0; j < size; j ++){
    mu += alpha/(alpha + std::pow(j, 1 - sigma));
    if (j == points_grid[id]) {
      res(id) = mu;
      id += 1;
    }
  }
  return(res);
}

// [[Rcpp::export]]
NumericVector polyseries_var_grid(NumericVector points_grid, int size, double alpha, double sigma) {
  int tot = points_grid.size();
  int id = 0;
  double v = 0.0;
  double term = 0.0;
  NumericVector res(tot);
  double mu = 0;
  for(int j = 0; j < size; j ++){
    term = alpha/(alpha + std::pow(j, 1 - sigma));
    v += term * (1 - term);
    if (j == points_grid[id]) {
      res(id) = v;
      id += 1;
    }
  }
  return(res);
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


