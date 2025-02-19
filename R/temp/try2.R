library(HubbellGLM)
X <- model.matrix(y ~ . - n, data = BarroColorado)
y <- BarroColorado$y
n <- BarroColorado$n

library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/polyseries.cpp")

polyseries_mean <- Vectorize(polyseries_mean_cpp, vectorize.args =c("size", "alpha"))
polyseries_var <- Vectorize(polyseries_var_cpp, vectorize.args = c("size", "alpha"))

g_inv_hub <- function(mu_target, size, sigma) {
  if (mu_target <= 1) {
    return(1e-7)
  } else if (mu_target >= size) {
    return(1e8)
  } else {
    uniroot(function(x) polyseries_mean(size, x, sigma) - mu_target,
            c(1e-10, 1e8), tol = 1e-8)$root
  }
}
g_inv_hub <- Vectorize(g_inv_hub, vectorize.args = c("mu_target", "size"))

mu <- 10
size <- 1000
sigma <- 0.1
g_inv_hub(mu_target = mu, size =  size, sigma = sigma)
uniroot(function(x) polyseries_mean(size, x, sigma) - mu, c(1e-10, 1e8), )

sigma * mu/(size-1)^sigma

Hubbell_NR(X, y, n, sigma = 0)$beta[, 1]
cbind(coef(fit_sim), Hubbell_NR(X, y, n, sigma = 0)$beta[, 1])

sigma <- 0.1
fit_sim <- HubbellGLMpoly(cbind(n, y) ~ ., data = BarroColorado,
                          family =hubbell.poly(sigma=0))
summary(fit_sim)
coef(fit_sim)
alphas <- g_inv_hub(mu_target = y, size = n, sigma = sigma)
coef(lm(log(alphas) ~ X - 1))

Hubbell_NR <- function(X, y, n, sigma = 0,
                       tol = 1e-16,
                       beta_start = NULL,
                       maxiter = 10000) {

  loglik <- numeric(maxiter)

  # Initialization (If null, implicitely initialized at beta=0)

  # Initialization
  mu <- y
  alpha <- g_inv_hub(mu, n, sigma = sigma)
  eta <- log(alpha)
  v <- alpha * (digamma(alpha + n) - digamma(alpha)) + alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
  w <- (1 / v) * polyseries_var(size =  n, alpha = alpha, sigma = sigma) ^ 2
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
    mu <- polyseries_mean(size = n, alpha = exp_eta, sigma = sigma)
    # Calculate the variance of the Antoniak distribution (does not depend on link, only on mu)
    alpha <- g_inv_hub(mu_target = mu, size = n, sigma = 0)
    v <- alpha * (digamma(alpha + n) - digamma(alpha)) +
      alpha^2 * (trigamma(alpha + n) - trigamma(alpha))
    # Calculate glm weights
    w <- (1 / v) * polyseries_var(size =  n, alpha = exp_eta, sigma = sigma) ^ 2
    # Calculate linearized response
    z <- eta + (y - mu) / sqrt(w * v)
    # Update loglikelihood
    loglik[t] <- sum(y * log(alpha) - lgamma(alpha + n) + lgamma(alpha))
    difference <- abs(loglik[t] - loglik[t - 1])
    #print(difference)
    if (difference < tol) {
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

res_tommi <- fisher_NR(X, y, n, maxiter = 10)
res_tommi$beta
res_tommi$Convergence[, 2]

res_ale <- Hubbell_NR(X, y, n, sigma = 0.1, maxiter = 10)
res_ale$Convergence
cbind("tommi" = res_tommi$beta[, 1], "ale" = res_ale$beta[, 1])


c("tommi_loglik" = tail(res_ale$Convergence, 1)[, 2],
  "ale_loglik" = tail(res_ale$Convergence, 1)[, 2])
res_tommi$Convergence

plot(X %*% res_ale$beta[, 1], X %*% res_tommi$beta[, 1])
abline(a = 0, b = 1)

beta <- res_tommi$beta[, 1]
exp_eta <- exp(X %*% beta)
mu <- polyseries_mean(n, alpha = exp_eta, 0)
alpha <- g_inv_hub(mu_target = mu, n, sigma = 0)
plot(exp_eta, alpha)
abline(a = 0, b = 1)
sum(y * log(alpha) - lgamma(alpha + n) + lgamma(alpha))



sourceCpp("src/test.cpp")
mu <- 87
sigma <- 0
n <- 10000

# Solve for alpha using Brent's method
a <- solve_alpha(mu, sigma, n)
b <- g_inv_hub(mu_target = mu, size =  n, sigma = sigma)
a
b

polyseries_mean(n, a, sigma) - mu
polyseries_mean(n, b, sigma) - mu

system.time( replicate(10000, solve_alpha_brent(mu, sigma, n)))
system.time( replicate(10000, g_inv_hub(mu_target = mu, size =  n, sigma = sigma)))

solve_alpha(100)

nsamples <- 500
mu <- runif(nsamples, min = 20, max = 1e3)
size <- sample(1e3:1e5, nsamples)
tt <- Sys.time()
res1 <- inv_polyseries(mu = mu, size = size, sigma = 0)
Sys.time() - tt
tt <- Sys.time()
res2 <- g_inv_hub(mu, size, sigma)
Sys.time() - tt

plot(res1, res2)
cbind(res1, res2)





