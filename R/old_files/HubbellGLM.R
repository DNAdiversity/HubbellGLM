#' Hubbell regression
#' This function mimics the GLM structure. We do so to avoid conflicts with the
#' current glm implementation (if it gets updated)
#'
#' @param formula Same as \code{stats::glm}
#' @param control Same as \code{stats::glm}
#' @param family a member of the hubbell family for \code{glm.fit.hubbell}.
#' @param data Same as \code{stats::glm}
#' @param weights Same as \code{stats::glm}
#' @param subset Same as \code{stats::glm}
#' @param na.action Same as \code{stats::glm}
#' @param start Same as \code{stats::glm}
#' @param etastart Same as \code{stats::glm}
#' @param mustart Same as \code{stats::glm}
#' @param offset Same as \code{stats::glm}
#' @param model Same as \code{stats::glm}
#' @param x Same as \code{stats::glm}
#' @param y Same as \code{stats::glm}
#' @param singular.ok Same as \code{stats::glm}
#' @param contrasts Same as \code{stats::glm}
#' @param ...Same as \code{stats::glm}
#'
#' @return An object of class "glm.hubell", which inherits from the class "lm"
#' @export
#'
HubbellGLM <- function(formula,
                       family = hubbell,
                       data,
                       weights,
                       subset,
                       na.action,
                       start = NULL,
                       etastart,
                       mustart,
                       offset,
                       control = list(...),
                       model = TRUE,
                       x = FALSE,
                       y = TRUE,
                       singular.ok = TRUE,
                       contrasts = NULL, ...) {
  # Force the method to be the glm for hubbell
  method <- "glm.fit.hubbell"

  # Here, the default glm code
  cal <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data)) {
    data <- environment(formula)
  }
  mf <- match.call(expand.dots = FALSE)
  m <- match(c(
    "formula", "data", "subset", "weights", "na.action",
    "etastart", "mustart", "offset"
  ), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame")) {
    return(mf)
  }
  if (!is.character(method) && !is.function(method)) {
    stop("invalid 'method' argument")
  }
  if (identical(method, "glm.fit")) {
    control <- do.call("glm.control", control)
  }
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) {
      names(Y) <- nm
    }
  }
  X <- if (!is.empty.model(mt)) {
    model.matrix(mt, mf, contrasts)
  } else {
    matrix(, NROW(Y), 0L)
  }
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) {
    stop("'weights' must be a numeric vector")
  }
  if (!is.null(weights) && any(weights < 0)) {
    stop("negative weights not allowed")
  }
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) {
      stop(gettextf(
        "number of offsets is %d should equal %d (number of observations)",
        length(offset), NROW(Y)
      ), domain = NA)
    }
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  fit <- eval(call(if (is.function(method)) "method" else method,
    x = X, y = Y, weights = weights, start = start, etastart = etastart,
    mustart = mustart, offset = offset, family = family,
    control = control, intercept = attr(mt, "intercept") >
      0L, singular.ok = singular.ok
  ))
  if (length(offset) && attr(mt, "intercept") > 0L) {
    fit2 <- eval(call(if (is.function(method)) "method" else method,
      x = X[, "(Intercept)", drop = FALSE], y = Y, mustart = fit$fitted.values,
      weights = weights, offset = offset, family = family,
      control = control, intercept = TRUE
    ))
    if (!fit2$converged) {
      warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
    }
    fit$null.deviance <- fit2$deviance
  }
  if (model) {
    fit$model <- mf
  }
  fit$na.action <- attr(mf, "na.action")
  if (x) {
    fit$x <- X
  }
  if (!y) {
    fit$y <- NULL
  }

  # Return an object of class
  structure(c(fit, list(
    call = cal, formula = formula, terms = mt,
    data = data, offset = offset, control = control, method = method,
    contrasts = attr(X, "contrasts"),
    xlevels = .getXlevels(mt, mf)
  )), class = c(fit$class, c("HubbellGLM", "lm")))
}


#' @export
family.HubbellGLM <- function(object, ...) {
  object$family
}

#' @export
logLik.HubbellGLM <- function(object, ...) {
  if (!missing(...)) {
    warning("extra arguments discarded")
  }
  fam <- family.HubbellGLM(object)$family
  dispersion <- family.HubbellGLM(object)$dispersion
  p <- object$rank
  if (!is.null(dispersion)) {
    if (is.na(dispersion)) {
      p <- p + 1
    }
  } else if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) {
    p <- p + 1
  }
  val <- p - object$aic / 2
  attr(val, "nobs") <- sum(!is.na(object$residuals))
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}
