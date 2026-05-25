#------------------------------ Hubbell Newton rapson
g_inv_hub <- function(mu_target, size, sigma) {
  if (mu_target <= 1) {
    return(1e-7)
  } else if (mu_target >= size) {
    return(1e8)
  } else {
    uniroot(function(x) HubbellGLM:::polyseries_mean(size, x, sigma) - mu_target,
            c(1e-10, 1e8), tol = 1e-8)$root
  }
}
g_inv_hub <- Vectorize(g_inv_hub, vectorize.args = c("mu_target", "size"))


Hubbell_NR <- function(X, y, n, sigma = 0,
                       tol = 1e-12,
                       beta_start = NULL,
                       maxiter = 10000) {

  loglik <- numeric(maxiter)

  # Initialization (If null, implicitely initialized at beta=0)

  # Initialization
  mu <- y
  alpha <- g_inv_hub(mu, n, sigma = sigma)
  eta <- log(alpha)
  v <- alpha * (digamma(alpha + n) - digamma(alpha)) + alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
  w <- (1 / v) * HubbellGLM:::polyseries_var(size =  n, alpha = alpha, sigma = sigma) ^ 2
  z <- eta + (y - mu) / sqrt(w * v)

  # First value of the likelihood
  loglik[1] <- sum(y * eta - lgamma(alpha + n) + lgamma(alpha))

  # Iterative procedure
  for (t in 2:maxiter) {
    # Find coefficients
    beta <- solve(qr(crossprod(X * w, X)), crossprod(X * w, z))
    eta <- c(X %*% beta)
    exp_eta <- exp(eta)
    # Calculate the mean
    mu <- HubbellGLM:::polyseries_mean(size = n, alpha = exp_eta, sigma = sigma)
    # Calculate the variance of the Antoniak distribution (does not depend on link, only on mu)
    alpha <- g_inv_hub(mu_target = mu, size = n, sigma = 0)
    v <- alpha * (digamma(alpha + n) - digamma(alpha)) +
      alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
    # Calculate glm weights
    w <- (1 / v) * HubbellGLM:::polyseries_var(size =  n, alpha = exp_eta, sigma = sigma) ^ 2
    # Calculate linearized response
    z <- eta + (y - mu) / sqrt(w * v)
    # Update loglikelihood
    loglik[t] <- sum(y * log(alpha) - lgamma(alpha + n) + lgamma(alpha))
    difference <- abs(loglik[t] - loglik[t - 1])

    if (difference < tol) {
      # Calculate fisher info and pvalues
      Fisher_info <- solve(t(X) %*% diag(w) %*% X)
      t_vals <- beta/sqrt(diag(Fisher_info))
      p_vals <- 2 * pnorm(-abs(t_vals))
      # Output the results
      results <- list(beta = cbind("beta" = beta,
                                   "se"= sqrt(diag(Fisher_info)),
                                   "t_vals" = t_vals,
                                   "p_vals" = p_vals),
                      vcov = Fisher_info,
                      loglik = loglik[t],
                      Convergence = cbind(Iteration = (1:t) - 1, Loglikelihood = loglik[1:t]))
      return(results)
    }
  }
  stop("The algorithm has not reached convergence")
}

#------------------------------ Fisher method (tommi's original function)
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
      # Calculate fisher info and pvalues
      Fisher_info <- solve(t(X) %*% diag(w) %*% X)
      t_vals <- beta/sqrt(diag(Fisher_info))
      p_vals <- 2 * pnorm(-abs(t_vals))
      # Output the results
      results <- list(beta = cbind(beta = beta,
                                   t_vals = t_vals,
                                   p_vals = p_vals),
                      vcov = Fisher_info,
                      loglik = loglik[t],
                      Convergence = cbind(Iteration = (1:t) - 1, Loglikelihood = loglik[1:t]))
      return(results)
    }
  }
  stop("The algorithm has not reached convergence")

}


#------------------------------ Simulation functions
simulate_data <- function(N, p, lambda, sigma = 0){
  X <- cbind(1, matrix(rnorm(N * (p - 1), sd = 10), nrow = N, ncol = p - 1))
  n <- rpois(N, lambda = lambda)
  beta <- c(2, rnorm(p-1, sd = 0.1))
  exp_eta <- exp(X %*% beta)
  mu <- HubbellGLM:::polyseries_mean(n, exp_eta, sigma  = sigma)
  # Find associated alphas
  alphas <- g_inv(mu, n)
  y <- rDPspecies(alpha = alphas, size = n)
  return(list(n = n, y = y,
              beta = beta, X = X,
              sigma= sigma,
              df = data.frame(n = n, y = y, X)))
}










