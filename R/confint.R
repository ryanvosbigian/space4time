
#' Confidence intervals for fitted model parameters for s4t_cjs
#'
#' @description
#' Computes Wald confidence intervals for one or more parameters in a fitted
#'     model object
#'
#' @param object A fitted model object of class `s4t_cjs`
#' @param parm A specification of what parameters are to be
#'    given confidence intervals. A character vector or names.
#' @param level The confidence level. Default is `0.95`
#' @param ... additional args. Not used
#' @return Gaussian-based confidence intervals for the fixed
#'      effect coffiecients based on the confidence level
#'      specified by `level`. Confidence intervals are on the
#'      logit scale.
#' @export
confint.s4t_cjs <- function(object, parm, level = 0.95, ...) {
  stopifnot(is(object,"s4t_cjs"))

  alpha <- 1 - level
  # tstar <- qt(1 - alpha / 2, df = object$n - object$p)
  tstar <- stats::qnorm(1 - alpha / 2)
  estimates <- object$estimated_parameters$estimate
  se <- object$estimated_parameters$std_error
  lower <- estimates - tstar * se
  upper <- estimates + tstar * se
  confints <- cbind(lower, upper)
  rownames(confints) <- object$estimated_parameters$parameter
  colnames(confints) <- c(paste(alpha / 2 * 100, "%"), paste((1 - alpha / 2) * 100, "%"))
  if (missing(parm)) {
    return(confints)
  } else {
    return(confints[row.names(confints) %in% parm, , drop = FALSE])
  }
}

