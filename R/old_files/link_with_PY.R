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

#------------------------------------------------ Pitman-yor link (polynomial)
#' Mean of the n. of species from a Pitman-Yor process
#'
#' @param alpha Precision parameter of the Pitman-Yor process
#' @param sigma Strength parameter of the Pitman-Yor process
#' @param n Number of units
#'
#' @export
mean_pitman_yor <- function(alpha, sigma, n){
  if(sigma == 0) {
    res <- mean_dirichlet_process(alpha, n)
  } else if (all(alpha > -sigma) & sigma > 0 & sigma < 1) {
    x <- log(alpha + n) + lgamma(alpha + sigma + n) - lgamma(alpha + sigma) -
      lgamma(alpha + n + 1) + lgamma(alpha + 1)
    res <- exp(x)/sigma - alpha/sigma
    #checks <- res > n | res < 1
    #res[checks] <- n[checks] - 1e-6
  } else {
    stop("Must have alpha > -sigma and sigma between (0, 1)")
  }
  return(res)
}


#' Derivative of the mean of the Dirichlet process distribution
#' @export
deriv_dirichlet_process <- function(alpha, n) {
  digamma(alpha + n) - digamma(alpha) + alpha * (trigamma(alpha + n) - trigamma(alpha))
}

#' Derivative of the mean of the Poisson-Dirichlet distribution
#' @export
deriv_pitman_yor <- function(alpha, sigma, n){
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
#' @export
inv_mean_dirichlet_process <- Vectorize(inv_mean_dirichlet_process, vectorize.args = c("mu_target", "size"))


# Find the associated alpha from the mean of the Poisson Dirichlet distribution
inv_mean_pitman_yor <- function(mu_target, size, sigma){
  if (mu_target <= 1 + 1e-4) {
    #sol <- .Machine$double.eps
    mu_target <- 1 + 1e-4
  } else if (mu_target >= size - 1e-4) {
    mu_target <- size - 1e-4
    #sol <- 1e18
  }
  if(sigma > 0 & sigma < 1) {
      # Solve case sigma in (0, 1),
      gamma <- uniroot(function(x) mean_pitman_yor(x - sigma, sigma, size) - mu_target,
                       c(.Machine$double.eps, 1e17), tol = 1e-6)$root
      alpha <- gamma - sigma
    } else if (sigma == 0){
      alpha <- inv_mean_dirichlet_process(mu_target, size)
    } else {
      stop("Must have 0 <= sigma < 1.")
    }
  return(alpha)
}
#' @export
inv_mean_pitman_yor <- Vectorize(inv_mean_pitman_yor, vectorize.args = c("mu_target", "size"))

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
         # Pitman-Yor mean
         py = {
           linkfun <- function(mu, size., sigma.) {
             log(inv_mean_pitman_yor(mu, size., sigma.) + sigma.)
           }
           # Inverse link function
           linkinv <- function(eta, size., sigma.) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             pmax(mean_pitman_yor(exp_eta - sigma., sigma., size.), .Machine$double.eps)
           }
           # Derivative of link function
           mu.eta <- function(eta, size., sigma.) {
             exp_eta <- pmax(exp(eta), .Machine$double.eps)
             exp_eta * deriv_pitman_yor(exp_eta - sigma., sigma., size.)
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


#' Hubbell family of genetalized linear models
#'
#' @param link  The link function to use. Available options are \code{dp}, \code{py} and \code{polyseries}
#' @param sigma Hyperparameter determining the polynomial growth.
#'              Must be 0 < sigma < 1 for \code{link = 'py'} and sigma < 1 for
#'              link  \code{link = 'polyseries'}. It is set to 0 when \code{link = 'dp'} (logarithmic growth).
#' @importFrom stats family
#' @export
hubbell <- function(link = "polyseries", sigma = 0) {

  # Call the link
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("dp", "py", "polyseries")
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
      stop("Must have sigma < 1 in link = c('py', 'polyseries')")
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
#'              Must be 0 < sigma < 1 for \code{link = 'py'} and sigma < 1 for
#'              link  \code{link = 'polyseries'}. It is set to 0 when \code{link = 'dp'} (logarithmic growth).
#' @importFrom stats family
#' @export
quasihubbell <- function(link = "polyseries", sigma = 0) {

  # Call the link
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  okLinks <- c("dp", "py", "polyseries")
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
      stop("Must have sigma < 1 in link = c('py', 'polyseries')")
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
