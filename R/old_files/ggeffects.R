




ggeffects_get_predictions_hubbell <- function (model, data_grid = NULL, terms = NULL, ci_level, type = NULL,
                                               typical = NULL, vcov = NULL, vcov_args = NULL, condition = NULL,
                                               interval = "confidence", bias_correction = FALSE, link_inverse = insight::link_inverse(model),
                                               model_info = NULL, verbose = TRUE, ...)
{
  se <- !is.null(ci_level) && !is.na(ci_level) && is.null(vcov)
  if (type == "simulate") {
    if (!is.null(ci_level) && !is.na(ci_level)) {
      ci <- (1 + ci_level)/2
    }
    else {
      ci <- 0.975
    }
    ggeffects:::.do_simulate(model, terms, ci, interval = interval, ...)
  }
  else {
    prdat <- suppressWarnings(predict.HubbellGLM(model, newdata = data_grid,
                                                 type = "link", se.fit = se, ...))
    ggeffects:::.generic_prediction_data(model, data_grid, link_inverse,
                             prediction_data = prdat, se, ci_level, typical, terms,
                             vcov, vcov_args, condition, interval)
  }
}



ggeffects_helper_hubbell <- function (model, terms, ci_level, type, typical, condition, back_transform,
          vcov, vcov_args, interval, bias_correction = FALSE, verbose = TRUE,
          ...)
{
  terms <- ggeffects:::.check_vars(terms, model)
  cleaned_terms <- ggeffects:::.clean_terms(terms)
  model_info <- ggeffects:::.get_model_info(model)
  if (inherits(model, "coxph") && type == "survival") {
    model_info$is_binomial <- TRUE
  }
  ggeffects:::.check_focal_for_random(model, terms, type, verbose)
  model_frame <- ggeffects:::.get_model_data(model)
  data_grid <- ggeffects:::.data_grid(model = model, model_frame = model_frame,
                          terms = terms, typical = typical, condition = condition,
                          show_pretty_message = verbose, verbose = verbose)
  original_model_frame <- model_frame
  original_terms <- terms
  terms <- cleaned_terms
  linv <- ggeffects:::.link_inverse(model, bias_correction = bias_correction,
                        ...)
  if (is.null(linv))
    linv <- function(x) x
  prediction_data <- ggeffects_get_predictions_hubbell(model, data_grid = data_grid,
                                     terms = original_terms, ci_level = ci_level, type = type,
                                     typical = typical, vcov = vcov, vcov_args = vcov_args,
                                     condition = condition, interval = interval, bias_correction = bias_correction,
                                     link_inverse = linv, model_info = model_info, verbose = verbose,
                                     ...)
  if (is.null(prediction_data)) {
    return(NULL)
  }
  attr(prediction_data, "continuous.group") <- attr(data_grid,
                                                    "continuous.group")
  if (inherits(model, "coxph") && type %in% c("survival", "cumulative_hazard")) {
    terms <- c("time", terms)
    cleaned_terms <- c("time", cleaned_terms)
  }
  if (inherits(model, "rqs") && !"tau" %in% cleaned_terms) {
    cleaned_terms <- c(cleaned_terms, "tau")
  }
  result <- ggeffects:::.post_processing_predictions(model = model, prediction_data = prediction_data,
                                         original_model_frame = original_model_frame, cleaned_terms = cleaned_terms)
  if (type == "simulate") {
    attributes(data_grid)$constant.values <- NULL
  }
  ggeffects:::.post_processing_labels_and_data(model = model, result = result,
                                   original_model_frame = original_model_frame, data_grid = data_grid,
                                   cleaned_terms = cleaned_terms, original_terms = original_terms,
                                   model_info = model_info, type = type, prediction.interval = attr(prediction_data,
                                                                                                    "prediction.interval", exact = TRUE), at_list = .data_grid(model = model,
                                                                                                                                                               model_frame = original_model_frame, terms = original_terms,
                                                                                                                                                               typical = typical, condition = condition, show_pretty_message = FALSE,
                                                                                                                                                               emmeans_only = TRUE, verbose = FALSE), condition = condition,
                                   ci_level = ci_level, back_transform = back_transform,
                                   vcov_args = .get_variance_covariance_matrix(model, vcov,
                                                                               vcov_args, skip_if_null = TRUE, verbose = FALSE),
                                   margin = "mean_reference", model_name = NULL, bias_correction = bias_correction,
                                   verbose = verbose)
}

#' @import ggeffects
#' @export
ggpredict_hubbell <- function (model, terms, ci_level = 0.95, type = "fixed", typical = "mean",
          condition = NULL, interval = "confidence", back_transform = TRUE,
          vcov = NULL, vcov_args = NULL, bias_correction = FALSE, verbose = TRUE,
          ...)
{
  type <- ggeffects:::.validate_type_argument(model, type)
  insight::formula_ok(model, verbose = verbose)
  interval <- insight::validate_argument(interval, c("confidence",
                                                     "prediction"))
  model.name <- deparse(substitute(model))
  bias_correction <- ggeffects:::.check_bias_correction(model, type = type,
                                            bias_correction = bias_correction, verbose = verbose)
  if (!missing(terms)) {
    terms <- ggeffects:::.reconstruct_focal_terms(terms, model = NULL)
  }
  model <- ggeffects:::.check_model_object(model)
  fun_args <- list(ci_level = ci_level, type = type, typical = typical,
                   condition = condition, back_transform = back_transform,
                   vcov = vcov, vcov_args = vcov_args, interval = interval,
                   bias_correction = bias_correction, verbose = verbose)
  if (inherits(model, "list") && !inherits(model, c("bamlss",
                                                    "maxLik"))) {
    result <- lapply(model, function(model_object) {
      full_args <- c(list(model = model_object, terms = terms),
                     fun_args, list(...))
      do.call(ggeffects_helper_hubbell, full_args)
    })
    class(result) <- c("ggalleffects", class(result))
  }
  else if (missing(terms) || is.null(terms)) {
    predictors <- insight::find_predictors(model, effects = "fixed",
                                           component = "conditional", flatten = TRUE, verbose = FALSE)
    result <- lapply(predictors, function(focal_term) {
      full_args <- c(list(model = model, terms = focal_term),
                     fun_args, list(...))
      tmp <- do.call(ggeffects_helper_hubbell, full_args)
      tmp$group <- focal_term
      tmp
    })
    names(result) <- predictors
    class(result) <- c("ggalleffects", class(result))
  }
  else {
    full_args <- c(list(model = model, terms = terms), fun_args,
                   list(...))
    result <- do.call(ggeffects_helper_hubbell, full_args)
  }
  if (!is.null(result)) {
    attr(result, "model.name") <- model.name
  }
  result
}



