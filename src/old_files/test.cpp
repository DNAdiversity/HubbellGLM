#include <Rcpp.h>
#include <boost/math/tools/roots.hpp>

using namespace Rcpp;
using namespace boost::math::tools;

// Function to evaluate f(alpha)
struct AlphaSolver {
  double mu;
  double sigma;
  int n;

  AlphaSolver(double mu_, double sigma_, int n_) : mu(mu_), sigma(sigma_), n(n_) {}

  // Function to compute f(alpha)
  double operator()(double alpha) const {
    double sum = 0.0;
    for (int j = 1; j <= n; j++) {
      double denom = alpha + std::pow(j - 1, 1 - sigma);
      sum += alpha / denom;
    }
    return sum - mu;
  }
};

// [[Rcpp::export]]
double solve_alpha(double mu, double sigma, int n, int max_iter = 1000) {
  double lower = 1e-16;  // Lower bound
  double upper = 1e16;   // Upper bound
  AlphaSolver solver(mu, sigma, n);
  // Tolerance and max iterations
  eps_tolerance<double> tol_func(40);  // Relative tolerance
  std::uintmax_t iters = max_iter;  // Track iteration count
  // Bisection method
  auto result = bisect(solver, lower, upper, tol_func, iters);
  // Extract the best root estimate (midpoint of the converged interval)
  double alpha = (result.first + result.second) / 2.0;
  return alpha;
}
