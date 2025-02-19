#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int DPspecies_cpp(double alpha, int size) {
  int k = 1;
  double u;
  for(int j = 1; j < size; j++){
    u = arma::randu();
    if(log(u) < log(alpha) -log(alpha + j)){
      k += 1;
    }
  }
  return k;
}

// [[Rcpp::export]]
arma::vec lastirling1(int n){

  arma::vec LogSk1(n+1); LogSk1.zeros();
  arma::vec LogSk(n+1);  LogSk.zeros();

  LogSk1(1) = 0;
  LogSk1(0) = -arma::datum::inf;
  LogSk(0)  = -arma::datum::inf;

  for(int i = 2; i <= n; i++){
    for(int j  = 1; j < i; j++){
      LogSk(j) = LogSk1(j) + std::log(i - 1 + std::exp(LogSk1(j-1) - LogSk1(j)));
    }
    LogSk(i)  = 0;
    LogSk1    = LogSk;
  }
  return(LogSk.rows(1,n));
}

// [[Rcpp::export]]
arma::mat lastirlings1(int n){

  arma::mat LogS(n+1,n+1);  LogS.fill(-arma::datum::inf);

  // Fill the starting values
  LogS(0,0) = 0;
  LogS(1,1) = 0;

  for(int i = 2; i <= n; i++){
    for(int j = 1; j < i; j++){
      LogS(i,j) = LogS(i-1,j) + std::log(i-1 + std::exp(LogS(i-1,j-1) - LogS(i-1,j)));
    }
    LogS(i,i)  = 0;
  }
  return(LogS.submat(1,1,n,n));
}










