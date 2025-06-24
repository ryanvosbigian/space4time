
#' Calculate variance-covariance matrix for a fitted model object
#'
#' @description
#' Calculate variance-covariance matrix for a `s4t_cjs_ml` model object.
#'
#'
#' @param object a `s4t_cjs_ml` object fitted using `fit_s4t_cjs_ml()`.
#' @param ... Other arguments. Not used (needed for generic consistency)
#'
#' @return The variance-covariance matrix of coefficients from a fitted
#'     `s4t_cjs_ml` object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim.dat <- sim_simple_s4t_ch(N = 2000)
#' m1 <- fit_s4t_cjs_ml(p_formula = ~ t,
#'                      theta_formula = ~ a1 * a2 * s * j,
#'                      ageclass_formula = ~ FL,
#'                      fixed_age = TRUE,
#'                      s4t_ch = sim.dat$s4t_ch)
#'
#' vcov(m1)
#' }
#'
vcov.s4t_cjs_ml <-  function(object, ...) {
  solve(object$res$hessian)
}



#' Calculate variance-covariance matrix for a fitted model object
#'
#' @description
#' Calculate variance-covariance matrix for a `s4t_ageclass_model` model object.
#'
#'
#' @param object a `s4t_ageclass_model` object fitted using `fit_ageclass()`.
#' @param ... Other arguments. Not used (needed for generic consistency)
#'
#' @return The variance-covariance matrix of coefficients from a fitted
#'     `s4t_ageclass_model` object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim.dat <- sim_simple_s4t_ch(N = 2000)
#' m1 <- fit_ageclass(ageclass_formula = ~ FL
#'                    s4t_ch = sim.dat$s4t_ch)
#'
#' vcov(m1)
#' }
#'
vcov.s4t_ageclass_model <-  function(object, ...) {
  solve(object$res$hessian)
}
