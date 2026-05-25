logEPPF_DP <- function(alpha, n, k) {
  # Loglikelihood
  loglik <- k * log(alpha) - lgamma(alpha + n) + lgamma(alpha)
  loglik
}

max_EPPF_DP <- function(n, k) {
  start <- 1 # Initialization of the maximization algorithm
  out <- nlminb(
    start = start,
    function(param) -logEPPF_DP(alpha = param, n = n, k = k),
    lower = 1e-10, upper = Inf
  )
  return(out$par)
}
max_EPPF_DP <- Vectorize(max_EPPF_DP, vectorize.args = c("n", "k"))

g_inv <- function(mu_target, size) {
  if (mu_target <= 1) {
    return(1e-7)
  } else if (mu_target >= size) {
    return(1e8)
  } else {
    uniroot(function(x) x * (digamma(x + size) - digamma(x)) - mu_target, c(1e-10, 1e8), tol = 1e-8)$root
  }
}
g_inv <- Vectorize(g_inv, vectorize.args = c("mu_target", "size"))

fisher_NR <- function(X, y, n,
                      tol = 1e-16,
                      beta_start = NULL,
                      maxiter = 10000) {
  loglik <- numeric(maxiter)
  # Initialization
  mu <- y
  alpha <- g_inv(mu, n)
  eta <- log(alpha)
  w <- mu + alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
  z <- eta + (y - mu) / w

  # First value of the likelihood
  loglik[1] <- sum(y * eta - lgamma(alpha + n) + lgamma(alpha))

  # Iterative procedure
  for (t in 2:maxiter) {
    beta <- solve(qr(crossprod(X * w, X)), crossprod(X * w, z))
    eta <- c(X %*% beta)
    alpha <- exp(eta)
    mu <- alpha * (digamma(alpha + n) - digamma(alpha))
    w <- mu + alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
    z <- eta + (y - mu) / w

    loglik[t] <- sum(y * eta - lgamma(alpha + n) + lgamma(alpha))
    if (abs(loglik[t] - loglik[t - 1]) < tol) {
      break
    }
  }
  #stop("The algorithm has not reached convergence")
  # Calculate fisher info and pvalues
  Fisher_info <- solve(t(X) %*% diag(w) %*% X)
  t_vals <- beta/sqrt(diag(Fisher_info))
  p_vals <- 2 * pnorm(-abs(t_vals))
  # Output the results
  results <- list(beta = cbind(beta = beta,
                               t_vals = t_vals,
                               p_vals = p_vals),
                  vcov = Fisher_info,
                  Convergence = cbind(Iteration = (1:t) - 1, Loglikelihood = loglik[1:t]))
  return(results)
}



