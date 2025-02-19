# Useful function that is present in the glm.fit...
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


# Function to fit the glm using the standard IRLS. Code is copied using the
# glm.fit function, adapting to include the parameter "size" in the code
glm.fit.hubbell.poly <- function(x, y,
                                 weights = rep.int(1, nobs),
                                 start = NULL,
                                 etastart = NULL,
                                 mustart = NULL,
                                 offset = rep.int(0, nobs),
                                 family = hubbell.poly(sigma = 0),
                                 control = list(),
                                 intercept = TRUE,
                                 singular.ok = TRUE) {
  control <- do.call("glm.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if (is.matrix(y)) {
    rownames(y)
  } else {
    names(y)
  }
  conv <- FALSE
  nobs <- NROW(y)
  nvars <- ncol(x)
  EMPTY <- nvars == 0
  if (is.null(weights)) {
    weights <- rep.int(1, nobs)
  }
  if (is.null(offset)) {
    offset <- rep.int(0, nobs)
  }
  variance <- family$variance
  linkinv <- family$linkinv
  if (!is.function(variance) || !is.function(linkinv)) {
    stop("'family' argument seems not to be a valid family object",
         call. = FALSE
    )
  }
  dev.resids <- family$dev.resids
  aic <- family$aic
  mu.eta <- family$mu.eta
  valideta <- family$valideta %||% function(eta) TRUE
  validmu <- family$validmu %||% function(mu) TRUE
  # Check the argument sigma for the polynomial series
  if(!is.null(family$sigma)){
    sigma <- family$sigma
  } else {
    sigma <- NULL
  }
  if (is.null(mustart)) {
    eval(family$initialize)
  } else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  if (EMPTY) {
    eta <- rep.int(0, nobs) + offset
    if (!valideta(eta)) {
      stop("invalid linear predictor values in empty model",
           call. = FALSE
      )
    }
    mu <- linkinv(eta, size, sigma)
    if (!validmu(mu)) {
      stop("invalid fitted means in empty model", call. = FALSE)
    }
    dev <- sum(dev.resids(y, mu, weights, size))
    w <- sqrt((weights * mu.eta(eta, size, sigma)^2) / variance(mu, size))
    residuals <- (y - mu) / mu.eta(eta, size, sigma)
    good <- rep_len(TRUE, length(residuals))
    boundary <- conv <- TRUE
    coef <- numeric()
    iter <- 0L
  } else {
    coefold <- NULL
    eta <- etastart %||% {
      if (!is.null(start)) {
        if (length(start) != nvars) {
          stop(
            gettextf(
              "length of 'start' should equal %d and correspond to initial coefs for %s",
              nvars, paste(deparse(xnames), collapse = ", ")
            ),
            domain = NA
          )
        } else {
          coefold <- start
          offset + as.vector(if (NCOL(x) == 1L) {
            x * start
          } else {
            x %*% start
          })
        }
      } else {
        family$linkfun(mustart, size, sigma)
      }
    }
    mu <- linkinv(eta, size, sigma)
    if (!(validmu(mu) && valideta(eta))) {
      stop("cannot find valid starting values: please specify some",
           call. = FALSE
      )
    }
    devold <- sum(dev.resids(y, mu, weights, size))
    boundary <- conv <- FALSE
    for (iter in 1L:control$maxit) {
      good <- weights > 0
      varmu <- variance(mu, size)[good]
      if (anyNA(varmu)) {
        stop("NAs in V(mu)")
      }
      if (any(varmu == 0)) {
        stop("0s in V(mu)")
      }
      mu.eta.val <- mu.eta(eta, size, sigma)
      if (any(is.na(mu.eta.val[good]))) {
        stop("NAs in d(mu)/d(eta)")
      }
      good <- (weights > 0) & (mu.eta.val != 0)
      if (all(!good)) {
        conv <- FALSE
        warning(gettextf(
          "no observations informative at iteration %d",
          iter
        ), domain = NA)
        break
      }
      z <- (eta - offset)[good] + (y - mu)[good] / mu.eta.val[good]
      w <- sqrt((weights[good] * mu.eta.val[good]^2) / variance(mu, size)[good])
      fit <- .Call(stats:::C_Cdqrls, x[good, , drop = FALSE] *
                     w, z * w, min(1e-07, control$epsilon / 1000), check = FALSE)
      if (any(!is.finite(fit$coefficients))) {
        conv <- FALSE
        warning(gettextf(
          "non-finite coefficients at iteration %d",
          iter
        ), domain = NA)
        break
      }
      if (nobs < fit$rank) {
        stop(sprintf(
          ngettext(
            nobs, "X matrix has rank %d, but only %d observation",
            "X matrix has rank %d, but only %d observations"
          ),
          fit$rank, nobs
        ), domain = NA)
      }
      if (!singular.ok && fit$rank < nvars) {
        stop("singular fit encountered")
      }
      start[fit$pivot] <- fit$coefficients
      eta <- drop(x %*% start)
      mu <- linkinv(eta <- eta + offset, size, sigma)
      dev <- sum(dev.resids(y, mu, weights, size))
      if (control$trace) {
        cat("Deviance = ", dev, " Iterations - ", iter,
            "\n",
            sep = ""
        )
      }
      boundary <- FALSE
      if (!is.finite(dev)) {
        if (is.null(coefold)) {
          stop("no valid set of coefficients has been found: please supply starting values",
               call. = FALSE
          )
        }
        warning("step size truncated due to divergence",
                call. = FALSE
        )
        ii <- 1
        while (!is.finite(dev)) {
          if (ii > control$maxit) {
            stop("inner loop 1; cannot correct step size",
                 call. = FALSE
            )
          }
          ii <- ii + 1
          start <- (start + coefold) / 2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset, size, sigma)
          dev <- sum(dev.resids(y, mu, weights, size))
        }
        boundary <- TRUE
        if (control$trace) {
          cat("Step halved: new deviance = ", dev, "\n",
              sep = ""
          )
        }
      }
      if (!(valideta(eta) && validmu(mu))) {
        if (is.null(coefold)) {
          stop("no valid set of coefficients has been found: please supply starting values",
               call. = FALSE
          )
        }
        warning("step size truncated: out of bounds",
                call. = FALSE
        )
        ii <- 1
        while (!(valideta(eta) && validmu(mu))) {
          if (ii > control$maxit) {
            stop("inner loop 2; cannot correct step size",
                 call. = FALSE
            )
          }
          ii <- ii + 1
          start <- (start + coefold) / 2
          eta <- drop(x %*% start)
          mu <- linkinv(eta <- eta + offset, size, sigma)
        }
        boundary <- TRUE
        dev <- sum(dev.resids(y, mu, weights, size))
        if (control$trace) {
          cat("Step halved: new deviance = ", dev, "\n",
              sep = ""
          )
        }
      }
      if (abs(dev - devold) / (0.1 + abs(dev)) < control$epsilon) {
        conv <- TRUE
        coef <- start
        break
      } else {
        devold <- dev
        coef <- coefold <- start
      }
    }
    if (!conv) {
      warning("glm.fit: algorithm did not converge", call. = FALSE)
    }
    if (boundary) {
      warning("glm.fit: algorithm stopped at boundary value",
              call. = FALSE
      )
    }
    eps <- 10 * .Machine$double.eps
    if (family$family == "binomial") {
      if (any(mu > 1 - eps) || any(mu < eps)) {
        warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                call. = FALSE
        )
      }
    }
    if (family$family == "poisson") {
      if (any(mu < eps)) {
        warning("glm.fit: fitted rates numerically 0 occurred",
                call. = FALSE
        )
      }
    }
    if (fit$rank < nvars) {
      coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
    }
    xxnames <- xnames[fit$pivot]
    residuals <- (y - mu) / mu.eta(eta, size, sigma)
    fit$qr <- as.matrix(fit$qr)
    nr <- min(sum(good), nvars)
    if (nr < nvars) {
      Rmat <- diag(nvars)
      Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
    } else {
      Rmat <- fit$qr[1L:nvars, 1L:nvars]
    }
    Rmat <- as.matrix(Rmat)
    Rmat[row(Rmat) > col(Rmat)] <- 0
    names(coef) <- xnames
    colnames(fit$qr) <- xxnames
    dimnames(Rmat) <- list(xxnames, xxnames)
  }
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  wt <- rep.int(0, nobs)
  wt[good] <- w^2
  names(wt) <- ynames
  names(weights) <- ynames
  names(y) <- ynames
  if (!EMPTY) {
    names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int(
      "",
      sum(good) - fit$rank
    ))
  }
  wtdmu <- if (intercept) {
    sum(weights * y) / sum(weights)
  } else {
    linkinv(offset, size, sigma)
  }
  nulldev <- sum(dev.resids(y, wtdmu, weights, size))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept)
  rank <- if (EMPTY) {
    0
  } else {
    fit$rank
  }
  resdf <- n.ok - rank
  aic.model <- aic(y, n, mu, weights, dev, size) + 2 * rank
  list(
    coefficients = coef, residuals = residuals, fitted.values = mu,
    effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
    rank = rank, qr = if (!EMPTY) {
      structure(fit[c(
        "qr", "rank",
        "qraux", "pivot", "tol"
      )], class = "qr")
    }, family = family,
    linear.predictors = eta, deviance = dev, aic = aic.model,
    null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
    df.residual = resdf, df.null = nulldf, y = y, size = size, sigma = sigma,
    converged = conv,
    boundary = boundary, name_size = name_size, name_y = name_y
  )
}
