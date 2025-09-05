
#' Efficient approximate leave-one-out cross-validation (LOO)
#'
#' Computes the PSIS-LOO CV efficient approximate leave-one-out (LOO)
#'     cross-validation for Bayesian models using Pareto smoothed importance sampling (PSIS).
#'     See documentation for `loo::loo()`
#'
#' @export
#'
#' @importFrom loo loo
#' @importFrom loo extract_log_lik
#'
#' @param x a `s4t_cjs_rstan` object
#' @param ... passed to loo::loo
#' @returns a `ADD DETAILS` object
loo.s4t_cjs_rstan <- function(x, ...) {

  log_lik_1 <- loo::extract_log_lik(x$res,parameter_name = "log_lik", merge_chains = FALSE)
  # loo::waic(log_lik_1)
  loo1 <- loo::loo(log_lik_1, ...)
  return(loo1)
}

#' Widely applicable information criterion (WAIC)
#'
#' add description
#' @export
#'
#' @importFrom loo waic
#' @importFrom loo extract_log_lik
#'
#' @param x a `s4t_cjs_rstan` object
#' @param ... not used
#' @returns A named list of class ....
waic.s4t_cjs_rstan <- function(x, ...) {

  log_lik_1 <- loo::extract_log_lik(x$res,parameter_name = "log_lik", merge_chains = FALSE)
  loo::waic(log_lik_1)
}
