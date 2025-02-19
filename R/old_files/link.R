# Link function for the Hubbell
#' @export
make.link.hubbell <- function(link) {
  switch(link,
    logseries = {
      linkfun <- function(mu, size. = size) {
        g_inv <- function(mu_target, s) {
          if (mu_target <= 1) {
            sol <- 1e-8
          } else if (mu_target >= s) {
            sol <- 1e8
          } else {
            sol <- uniroot(function(x) x * (digamma(x + s) - digamma(x)) - mu_target, c(1e-12, 1e12))$root
          }
          log(sol)
        }
        g_inv <- Vectorize(g_inv, vectorize.args = c("mu_target", "s"))
        g_inv(mu, size.)
      }

      linkinv <- function(eta, size. = size) {
        exp_eta <- pmax(exp(eta), .Machine$double.eps)
        pmax(exp_eta * (digamma(exp_eta + size.) - digamma(exp_eta)), .Machine$double.eps)
      }

      mu.eta <- function(eta, size. = size) {
        exp_eta <- pmax(exp(eta), .Machine$double.eps)
        pmax(exp_eta * (digamma(exp_eta + size.) - digamma(exp_eta)) +
          exp_eta^2 * (trigamma(exp_eta + size.) - trigamma(exp_eta)), .Machine$double.eps)
      }

      valideta <- function(eta, size. = size) TRUE
    },
    stop(gettextf("%s link not recognised", sQuote(link)),
      domain = NA
    )
  )
  environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(valideta) <- asNamespace("stats")
  structure(list(
    linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = "hubbell"
  ), class = "link-glm")
}

# Hubbell family
#' @export
hubbell <- function(link = "logseries") {
  family <- "hubbell"
  # Make the accumulation curve link for the Antoniak glm
  stats <- make.link.hubbell(link)
  linktemp <- "logseries"

  # Variance function
  variance <- function(mu, size. = size) {
    exp_eta <- exp(stats$linkfun(mu, size.))
    mu + exp_eta^2 * (trigamma(exp_eta + size.) - trigamma(exp_eta))
  }

  # Valid values for mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu >= 1)

  # residual deviance
  dev.resids <- function(y, mu, wt, size. = size) {
    alpha_y <- exp(stats$linkfun(y, size.))
    alpha_mu <- exp(stats$linkfun(mu, size.))
    2 * (wt * (y * log(alpha_y / alpha_mu) - lgamma(alpha_y + size.) + lgamma(alpha_y) +
      lgamma(alpha_mu + size.) - lgamma(alpha_mu)))
  }

  # aic
  aic <- function(y, n, mu, wt, dev, size. = size) {
    alpha <- exp(stats$linkfun(mu, size.))
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
    dispersion = 1
  ), class = "family")
}


#' @export
make.link.hubbell.poly <- function(link) {
  switch(link,
         polyseries = {
           linkfun <- function(mu, size. = size, sigma. = sigma) {
             g_inv <- function(mu_target, s, sigma) {
               if (mu_target <= 1) {
                 sol <- 1e-8
               } else if (mu_target >= s) {
                 sol <- 1e8
               } else {
                 sol <- uniroot(function(x) polyseries_mean(s, x, sigma) - mu_target, c(1e-12, 1e12))$root
               }
               log(sol)
             }
             g_inv <- Vectorize(g_inv, vectorize.args = c("mu_target", "s"))
             g_inv(mu, size., sigma.)
           }

           # Inverse link function
           linkinv <- function(eta, size. = size, sigma. = sigma) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             pmax(polyseries_mean(size., exp_eta, sigma.), .Machine$double.eps)
           }

           # Derivative of link function
           mu.eta <- function(eta, size. = size, sigma. = sigma) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             pmax(polyseries_var(size = size., alpha = exp_eta, sigma = sigma.), .Machine$double.eps)
           }

           valideta <- function(eta, size. = size) TRUE
         },
         stop(gettextf("%s link not recognised", sQuote(link)), domain = NA)
  )
  environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(valideta) <- asNamespace("stats")
  structure(list(
    linkfun = linkfun,
    linkinv = linkinv,
    mu.eta = mu.eta,
    valideta = valideta,
    name = "hubbell"
  ), class = "link-glm")
}

# Exponential family
#' @export
hubbell.poly <- function(link = "polyseries", sigma = 0) {
  family <- "hubbell.poly"
  # Make the accumulation curve link for the Hubbell glm
  stats <- make.link.hubbell.poly(link)
  linktemp <- "polyseries"

  # Variance function
  variance <- function(mu, size. = size) {
    exp_eta <- exp(stats$linkfun(mu, size., sigma. = 0))
    mu + exp_eta^2 * (trigamma(exp_eta + size.) - trigamma(exp_eta))
  }

  # Valid values for mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu >= 1)

  # residual deviance
  dev.resids <- function(y, mu, wt, size. = size) {
    alpha_y <- exp(stats$linkfun(y, size., sigma. = 0))
    alpha_mu <- exp(stats$linkfun(mu, size., sigma. = 0))
    2 * (wt * (y * log(alpha_y / alpha_mu) - lgamma(alpha_y + size.) + lgamma(alpha_y) +
                 lgamma(alpha_mu + size.) - lgamma(alpha_mu)))
  }

  # aic
  aic <- function(y, n, mu, wt, dev, size. = size) {
    alpha <- exp(stats$linkfun(mu, size., sigma. = 0))
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
    dispersion = 1, sigma = sigma), class = "family")
}


