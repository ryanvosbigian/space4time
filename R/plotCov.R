

#' Effect plot for s4t_cjs_rstan object
#'
#' add description
#'
#'
#' @param mod a `s4t_cjs_rstan` object
#' @param param a `character` of the parameter
#' @param ... not used at the moment
#' @returns a ggplot2 type figure of the transition probabilities
plotCov <- function(mod, param, ...) {
  if (!is(mod,"s4t_cjs_rstan")) stop("mod must be class s4t_cjs_ml or s4t_cjs_rstan")

  if (param %in% c("a1","a2","s","t","j","k","b","g")) stop("use plotTheta or plotSurvival for this param")

  # param # check that it is a character of length 1



  # if (mod$call$theta_formula)

  indices_theta_original <- mod$original_units$indices_theta_original
  indices_theta_incCov <- indices_theta_original[,1:(ncol(indices_theta_original) - 8)] %>%
    dplyr::mutate(a1 = factor(a1),
                  a2 = factor(a2),
                  s = factor(s),
                  t = factor(t),
                  j = factor(j),
                  k = factor(k),
                  b = factor(b),
                  g = factor(g))




  mod_mat = model_mat_info(mod$call$theta_formula,df = indices_theta_incCov)$mod_mat

  # 1. determine if param is a factor (being used as a factor), use mod_mat to decide
  # 2. if factor: loop through and manipulate mod_mat so that it turns on each level at a time.
  #     also need to do something about interactions.
  # 3. if numeric, loop through and manipulate mod_mat so that it changes the value one step at a time.
  #     will need to deal with interactions.
  # 4. make a figure. x-axis will be the covariate/param. The y-axis will be cohort_transition/surv.
  #     Add facet_wrap for factors and other variables?

  return(NULL)

}

