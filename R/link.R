# Polynomial series for hubbell regression

#' Function to invert the polyseries
#' @param mu_target Target mean
#' @param size Number of units in the sequence
#' @param sigma Parameter for the growth. sigma = 0 is the logarithmic growth
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib HubbellGLM
#' @export
inv_polyseries <- function(mu_target, size, sigma) {
  if (mu_target <= 1 + 1e-4) {
    #sol <- .Machine$double.eps
    mu_target <- 1 + 1e-4
  } else if (mu_target >= size - 1e-4) {
    mu_target <- size - 1e-4
    #sol <- 1e18
  }
  # Invert the function
  lb <- (mu_target - 1) / (size - mu_target)
  ub <- (size - 1) ^ (1 - sigma) * (mu_target - 1) / (size - mu_target) + 1e-4
  sol <- uniroot(function(x) HubbellGLM:::polyseries_mean(size, x, sigma) - mu_target,
                 c(lb, ub), tol = 1e-6)$root
  return(sol)
}

# Inverse function for the polynomial series
#' @export
inv_polyseries <- Vectorize(inv_polyseries, vectorize.args = c("mu_target", "size"))


#------------------------------------------------ Dirichlet process link (canonical)
#' Mean of the n. of species from Dirichlet process
#'
#' @param alpha parameter of the Dirichlet process
#' @param n Number of units
#'
#' @export
mean_dirichlet_process <- function(alpha, n) {
  res <- alpha * (digamma(alpha + n) - digamma(alpha))
  return(res)
}


#' Derivative of the mean of the Dirichlet process distribution
#'
#' @param alpha Concentration parameter (positive scalar or vector).
#' @param n Sample size (positive integer or vector).
#'
#' @export
deriv_dirichlet_process <- function(alpha, n) {
  digamma(alpha + n) - digamma(alpha) + alpha * (trigamma(alpha + n) - trigamma(alpha))
}

# Find the associated alpha from the mean of the Dirichlet process
inv_mean_dirichlet_process <- function(mu_target, size){
  if (mu_target <= 1 + 1e-4) {
    #sol <- .Machine$double.eps
    mu_target <- 1 + 1e-4
  } else if (mu_target >= size - 1e-4) {
    mu_target <- size - 1e-4
    #sol <- 1e18
  }
  lb <- (mu_target - 1) / (size - mu_target)
  ub <- (size - 1) * (mu_target - 1) / (size - mu_target) + 1e-4
  alpha <- uniroot(function(x) mean_dirichlet_process(x, size) - mu_target,
                   c(lb, ub), tol = 1e-6)$root
  return(alpha)
}
#' Invert the mean of the Dirichlet process distribution
#'
#' Given a target mean \eqn{\mu} and a community size \eqn{n}, finds the
#' concentration parameter \eqn{\alpha} such that
#' \eqn{E[Y^{(n)}; \alpha] = \mu} via \code{uniroot}.
#'
#' @param mu_target Target mean value (scalar or vector, must be in
#'   \eqn{(1, n)}).
#' @param size Community size (positive integer, scalar or vector).
#'
#' @return The concentration parameter \eqn{\alpha} (scalar or vector).
#'
#' @export
inv_mean_dirichlet_process <- Vectorize(inv_mean_dirichlet_process, vectorize.args = c("mu_target", "size"))

#' Link function for the Hubbell generalized linear model
#' @param link Name for the link function
#' @importFrom Rcpp sourceCpp
#' @useDynLib HubbellGLM
#' @export
make.link.hubbell <- function(link) {
  switch(link,
         # Polyseries link
         polyseries = {
           linkfun <- function(mu, size., sigma.) {
             log(inv_polyseries(mu, size., sigma.))
           }
           # Inverse link function
           linkinv <- function(eta, size., sigma.) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             pmax(HubbellGLM:::polyseries_mean(size., exp_eta, sigma.), .Machine$double.eps)
           }

           # Derivative of link function
           mu.eta <- function(eta, size., sigma.) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             pmax(HubbellGLM:::polyseries_var(size = size., alpha = exp_eta, sigma = sigma.),
                  .Machine$double.eps)
           }
           valideta <- function(eta, size.) TRUE
         },
         # Dirichlet process link (canonical)
         dp = {
           linkfun <- function(mu, size., sigma.) {
             log(inv_mean_dirichlet_process(mu, size.))
           }

           # Inverse link function
           linkinv <- function(eta, size., sigma.) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             pmax(mean_dirichlet_process(exp_eta, size.), .Machine$double.eps)
           }

           # Derivative of link function
           mu.eta <- function(eta, size., sigma.) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             exp_eta * deriv_dirichlet_process(exp_eta, size.)
           }
           valideta <- function(eta, size.) TRUE
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


#' Hubbell family of genetalized linear models
#'
#' @param link  The link function to use. Available options are \code{dp}, and \code{polyseries}
#' @param sigma Hyperparameter determining the polynomial growth. Must be \code{sigma < 1} for
#'              link  \code{link = 'polyseries'}. It is set to 0 when \code{link = 'dp'} (logarithmic growth).
#' @importFrom stats family
#' @export
hubbell <- function(link = "polyseries", sigma = 0) {

  # Call the link
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("dp", "polyseries")
  family <- "hubbell"
  if (linktemp %in% okLinks)
    stats <- make.link.hubbell(linktemp)
  else if (is.character(link)) {
    stats <- make.link.hubbell(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for %s family; available links are %s",
                    linktemp, family, paste(sQuote(okLinks), collapse = ", ")),
           domain = NA)
    }
  }
  if(linktemp %in% c("polyseries")) {
    if(sigma >= 1){
      stop("Must have sigma < 1 in link = c('polyseries')")
    }
    linktemp <- paste0("linktemp", " - sigma = ", sigma)
  }

  # Variance function
  variance <- function(mu, size.) {
    # Invert the function to find the concentration parameter
    alpha <- inv_mean_dirichlet_process(mu, size.) #inv_polyseries(mu, size., sigma = 0)
    mu + alpha^2 * (trigamma(alpha + size.) - trigamma(alpha)) #polyseries_var(size., alpha, 0) #
  }

  # Valid values for mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu >= 1)

  # residual deviance
  dev.resids <- function(y, mu, wt, size.) {
    alpha_y <- inv_mean_dirichlet_process(y, size.) #inv_polyseries(y, size., sigma = 0)
    alpha_mu <- inv_mean_dirichlet_process(mu, size.) #inv_polyseries(mu, size., sigma = 0)
    2 * (wt * (y * log(alpha_y / alpha_mu) - lgamma(alpha_y + size.) + lgamma(alpha_y) +
                 lgamma(alpha_mu + size.) - lgamma(alpha_mu)))
  }

  # aic
  aic <- function(y, n, mu, wt, dev, size.) {
    alpha <- inv_mean_dirichlet_process(mu, size.) #inv_polyseries(mu, size., sigma = 0)
    # Calculate the log signless stirling numbers of the first kind
    Stir <- lastirlings1(max(size.))
    Stir <- apply(cbind(size., y), 1, function(x) Stir[x[1], x[2]])
    -2 * sum((y * log(alpha) - lgamma(alpha + size.) + lgamma(alpha) + Stir) * wt)
  }

  # Initialization function
  initialize <- expression({
    #print(y)
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


#' Quasi-Hubbell family of generalized linear models
#'
#' @param link  The link function to use. Available options are \code{dp}, \code{py} and \code{polyseries}
#' @param sigma Hyperparameter determining the polynomial growth.
#'              Must be \code{sigma < 1} for link  \code{link = 'polyseries'}. It is set to 0 when \code{link = 'dp'} (logarithmic growth).
#' @importFrom stats family
#' @export
quasihubbell <- function(link = "polyseries", sigma = 0) {

  # Call the link
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("dp", "polyseries")
  family <- "hubbell"
  if (linktemp %in% okLinks)
    stats <- make.link.hubbell(linktemp)
  else if (is.character(link)) {
    stats <- make.link.hubbell(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for %s family; available links are %s",
                    linktemp, family, paste(sQuote(okLinks), collapse = ", ")),
           domain = NA)
    }
  }
  if(linktemp %in% c("py", "polyseries")) {
    if(sigma >= 1){
      stop("Must have sigma < 1 in link = c('polyseries')")
    }
    linktemp <- paste0("linktemp", " - sigma = ", sigma)
  }

  # Variance function
  variance <- function(mu, size.) {
    # Invert the function to find the concentration parameter
    #alpha <- inv_polyseries(mu, size., sigma = 0)
    alpha <- inv_mean_dirichlet_process(mu, size.)
    mu + alpha^2 * (trigamma(alpha + size.) - trigamma(alpha))
  }

  # Valid values for mu
  validmu <- function(mu) all(is.finite(mu)) && all(mu >= 1)

  # residual deviance
  dev.resids <- function(y, mu, wt, size.) {
    alpha_y <- inv_mean_dirichlet_process(y, size.)#inv_polyseries(y, size., sigma = 0)
    alpha_mu <- inv_mean_dirichlet_process(mu, size.) #inv_polyseries(mu, size., sigma = 0)
    2 * (wt * (y * log(alpha_y / alpha_mu) - lgamma(alpha_y + size.) + lgamma(alpha_y) +
                 lgamma(alpha_mu + size.) - lgamma(alpha_mu)))
  }

  # aic
  aic <- function(y, n, mu, wt, dev, size.) {
    alpha <- inv_mean_dirichlet_process(mu, size.) #inv_polyseries(mu, size., sigma = 0)
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
