


#' @export
quasihubbell <- function(link = "polyseries", sigma = 0) {
  family <- "hubbell"
  # Make the accumulation curve link for the Hubbell glm
  stats <- make.link.hubbell(link)
  linktemp <- paste0("polyseries", " - sigma = ", sigma)

  # Variance function
  variance <- function(mu, size.) {
    # Invert the function to find the correct polyseries
    alpha <- inv_polyseries(mu, size., sigma = 0)
    mu + alpha^2 * (trigamma(alpha + size.) - trigamma(alpha))
  }

  # Valid values for mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu >= 1)

  # residual deviance
  dev.resids <- function(y, mu, wt, size.) {
    alpha_y <- inv_polyseries(y, size., sigma = 0)
    alpha_mu <- inv_polyseries(mu, size., sigma = 0)
    2 * (wt * (y * log(alpha_y / alpha_mu) - lgamma(alpha_y + size.) + lgamma(alpha_y) +
                 lgamma(alpha_mu + size.) - lgamma(alpha_mu)))
  }

  # aic
  aic <- function(y, n, mu, wt, dev, size.) {
    alpha <- inv_polyseries(mu, size., sigma = 0)
    # Calculate the log signless stirling numbers of the first kind
    Stir <- lastirlings1(max(size.))
    Stir <- apply(cbind(size., y), 1, function(x) Stir[x[1], x[2]])
    -2 * sum((y * log(alpha) - lgamma(alpha + size.) + lgamma(alpha) + Stir) * wt)
  }

  # Initialization function
  initialize <- expression({
    if (ncol(y) != 2) stop("must specify cbind(n, y) as dependent variable for the 'hubbell' family")
    name_size <- colnames(y)[1]
    name_y <- colnames(y)[2]
    size <- y[, 1]
    mustart <- y[, 2]
    y <- y[, 2]
    if (any(y < 1 | y > size)) stop("values not satisfying 1 <= y <= size for the 'hubbell' family")
  })

  # Simulation function
  simfun <- function(object, nsim, size.) {
    wts <- object$prior.weights
    if (any(wts != 1)) {
      warning("ignoring prior weights")
    }
    ftd <- fitted(object)
    alpha <- exp(stats$linkfun(ftd))
    sapply(1:nsim, function(i) rDPspecies(alpha = alpha, size = size.))
  }

  # Return the structure for thr glm type
  structure(list(
    family = family, link = linktemp, linkfun = stats$linkfun,
    linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
    aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
    validmu = validmu, valideta = stats$valideta, simulate = simfun,
    dispersion = NA_real_, sigma = sigma), class = "family")
}


## Exponential family
## @export
# hubbell <- function(link = "polyseries", sigma = 0) {
#   family <- "hubbell"
#   # Make the accumulation curve link for the Hubbell glm
#   stats <- make.link.hubbell(link)
#   linktemp <- paste0("polyseries", " - sigma = ", sigma)
#
#   # Variance function
#   variance <- function(mu, size.) {
#     # Invert the function to find the correct polyseries
#     alpha <- inv_polyseries(mu, size., sigma = 0)
#     mu + alpha^2 * (trigamma(alpha + size.) - trigamma(alpha))
#   }
#
#   # Valid values for mu
#   validmu <- function(mu) all(is.finite(mu)) && all(mu >= 1)
#
#   # residual deviance
#   dev.resids <- function(y, mu, wt, size.) {
#     alpha_y <- inv_polyseries(y, size., sigma = 0)
#     alpha_mu <- inv_polyseries(mu, size., sigma = 0)
#     2 * (wt * (y * log(alpha_y / alpha_mu) - lgamma(alpha_y + size.) + lgamma(alpha_y) +
#                  lgamma(alpha_mu + size.) - lgamma(alpha_mu)))
#   }
#
#   # aic
#   aic <- function(y, n, mu, wt, dev, size.) {
#     alpha <- inv_polyseries(mu, size., sigma = 0)
#     # Calculate the log signless stirling numbers of the first kind
#     Stir <- lastirlings1(max(size.))
#     Stir <- apply(cbind(size., y), 1, function(x) Stir[x[1], x[2]])
#     -2 * sum((y * log(alpha) - lgamma(alpha + size.) + lgamma(alpha) + Stir) * wt)
#   }
#
#   # Initialization function
#   initialize <- expression({
#     if (ncol(y) != 2) stop("must specify cbind(n, y) as dependent variable for the 'hubbell' family")
#     name_size <- colnames(y)[1]
#     name_y <- colnames(y)[2]
#     size <- y[, 1]
#     mustart <- y[, 2]
#     y <- y[, 2]
#     if (any(y < 1 | y > size)) stop("values not satisfying 1 <= y <= size for the 'hubbell' family")
#   })
#
#   # Simulation function
#   simfun <- function(object, nsim, size.) {
#     wts <- object$prior.weights
#     if (any(wts != 1)) {
#       warning("ignoring prior weights")
#     }
#     ftd <- fitted(object)
#     alpha <- exp(stats$linkfun(ftd))
#     sapply(1:nsim, function(i) rDPspecies(alpha = alpha, size = size.))
#   }
#
#   # Return the structure for thr glm type
#   structure(list(
#     family = family, link = linktemp, linkfun = stats$linkfun,
#     linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
#     aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
#     validmu = validmu, valideta = stats$valideta, simulate = simfun,
#     dispersion = 1, sigma = sigma), class = "family")
# }

#---------------------------------------------------------------
# Old functions, not to use
#' Mean of the n. of species from Poisson-Dirichlet distribution
mean_poisson_dirichlet <- function(alpha, sigma, n){
  if(sigma == 0) {
    res <- mean_dirichlet_process(alpha, n)
  } else if (all(alpha > -sigma)) {
    x <- log(alpha + n) + lgamma(alpha + sigma + n) - lgamma(alpha + sigma) -
      lgamma(alpha + n + 1) + lgamma(alpha + 1)
    res <- exp(x)/sigma - alpha/sigma
    if(res > n | res < 1){
      # This indicates overflow.
      res <- n - 1e-6
    }
  } else {
    stop("Must have alpha > -sigma")
  }
  return(res)
}

deriv_poisson_dirichlet <- function(alpha, sigma, n){
  if(sigma == 0){
    dd <- deriv_dirichlet_process(alpha, n)
  } else {
    gamma <- alpha + sigma # The derivative is respect to the quanity gamma = alpha + sigma
    log_poch_r <- lgamma(gamma + n) - lgamma(gamma) -
      lgamma(alpha + n + 1) + lgamma(alpha + 1)
    # Return the derivative
    dd <- exp(log_poch_r)/sigma -
      1/sigma +
      exp(log_poch_r + log(alpha + n))/sigma * (digamma(gamma + n) -
                                                  digamma(gamma) -
                                                  digamma(alpha + n + 1) +
                                                  digamma(alpha + 1))
  }
  return(dd)
}


# Find the associated alpha from the mean of the Poisson Dirichlet distribution
inv_mean_poisson_dirichlet <- function(mu_target, size, sigma){
  if (mu_target <= 1 + 1e-6) {
    alpha <- -sigma + .Machine$double.eps
  } else if (mu_target >= size - 1e-6) {
    alpha <- 1e18
  } else {
    if(sigma > 0 & sigma < 1) {
      # Solve case sigma in (0, 1), using logarithms
      alpha <- uniroot(function(x) mean_poisson_dirichlet(x, sigma, size) - mu_target,
                       c(-sigma + .Machine$double.eps, 1e10), tol = 1e-6)$root
    }
  }
  return(alpha)
}
inv_mean_poisson_dirichlet <- Vectorize(inv_mean_poisson_dirichlet, vectorize.args = c("mu_target", "size"))




