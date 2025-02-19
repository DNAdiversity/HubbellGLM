# Function to predict
#' @export
predict.HubbellGLM <- function(
    object, newdata = NULL, type = c(
      "link", "response",
      "terms"
    ), se.fit = FALSE, dispersion = NULL, terms = NULL,
    na.action = na.pass, ...) {
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL
  if (!se.fit) {
    if (missing(newdata)) {
      pred <- switch(type,
        link = object$linear.predictors,
        response = object$fitted.values,
        terms = predict.lm(object,
          se.fit = se.fit, scale = 1, type = "terms",
          terms = terms
        )
      )
      if (!is.null(na.act)) {
        pred <- napredict(na.act, pred)
      }
    } else {
      pred <- predict.lm(object, newdata, se.fit,
        scale = 1,
        type = if (type == "link") {
          "response"
        } else {
          type
        }, terms = terms, na.action = na.action
      )
      switch(type,
        response = {
          if (object$family$family != "hubbell") {
            pred <- family(object)$linkinv(pred)
          } else {
            sigma <- family(object)$sigma
            pred <- family(object)$linkinv(pred, newdata[, object$name_size], sigma)
          }
        },
        link = ,
        terms =
        )
    }
  } else {
    if (inherits(object, "survreg")) {
      dispersion <- 1
    }
    if (is.null(dispersion) || dispersion == 0) {
      dispersion <- summary(object, dispersion = dispersion)$dispersion
    }
    residual.scale <- as.vector(sqrt(dispersion))
    pred <- predict.lm(object, newdata, se.fit,
      scale = residual.scale,
      type = if (type == "link") {
        "response"
      } else {
        type
      }, terms = terms, na.action = na.action
    )
    fit <- pred$fit
    se.fit <- pred$se.fit
    switch(type,
      response = {
        if (object$family$family != "hubbell") {
          se.fit <- se.fit * abs(family(object)$mu.eta(fit))
          fit <- family(object)$linkinv(fit)
        } else {
          sigma <- family(object)$sigma
          if (missing(newdata)) {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit, object$size, sigma))
            fit <- family(object)$linkinv(fit, object$size, sigma)
          } else {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit, newdata[, object$name_size], sigma))
            fit <- family(object)$linkinv(fit, newdata[, object$name_size], sigma)
          }
        }
      },
      link = ,
      terms =
      )
    if (missing(newdata) && !is.null(na.act)) {
      fit <- napredict(na.act, fit)
      se.fit <- napredict(na.act, se.fit)
    }
    pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
  }
  pred
}
