#' @export
residuals.HubbellGLM <- function(object, type = c("deviance", "pearson", "working", "response", "partial"), ...) {
  type <- match.arg(type)
  y <- object$y
  r <- object$residuals
  mu <- object$fitted.values
  wts <- object$prior.weights
  switch(type,
    deviance = ,
    pearson = ,
    response = if (is.null(y)) {
      mu.eta <- object$family$mu.eta
      eta <- object$linear.predictors
      y <- mu + r * mu.eta(eta)
    }
  )

  size <- object$size
  res <- switch(type,
    deviance = if (object$df.residual > 0) {
      d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, wts, size), 0))
      ifelse(y > mu, d.res, -d.res)
    } else {
      rep.int(0, length(mu))
    },
    pearson = (y - mu) * sqrt(wts) / sqrt(object$family$variance(mu, size)),
    working = r,
    response = y - mu,
    partial = r
  )

  if (!is.null(object$na.action)) {
    res <- naresid(object$na.action, res)
  }
  if (type == "partial") {
    res <- res + predict(object, type = "terms")
  }
  res
}
