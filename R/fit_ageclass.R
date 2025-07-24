

#' Fit age-class model using ordinal regression
#'
#' @description
#' Fit age-class ordinal regression model with logit link and
#'     flexible threshold parameters.
#'
#'
#' @export
#' @param age_formula a one- or two-sided `formula` describing the effect
#'     structure of the model.
#' @param s4t_ch a `s4t_ch` object
#' @returns a fitted `s4t_ageclass_model` object that contains the estimated
#'     parameters, AIC, negative log likelihood, and call.
#'
#' @examples
#' sim.dat <- sim_simple_s4t_ch()
#' m <- fit_ageclass(ageclass ~ FL + obs_time,s4t_ch = sim.dat$s4t_ch)
#' summary(m)
#'
fit_ageclass <- function(age_formula = ageclass ~ 1,
                         s4t_ch) {

  # validate_s4t_ch(s4t_ch)

  ### add checks

  new_max_a <- max(s4t_ch$ch$all_aux[,"ageclass"],na.rm = TRUE)

  ageclass_data <- ageclass_call(age_formula=age_formula,
                                 obs_aux = stats::na.omit(s4t_ch$ch$all_aux)
                                 # max_a = s4t_ch$ch_info$max_a
  )


  inits <- ageclass_data$age_par

  lower <- ifelse(grepl("a_alpha_",names(inits)) &
                    !grepl("^a_alpha_1$",names(inits)),
                  0,-Inf)

  upper <- rep(Inf,length(inits))

  names(lower) <- names(upper) <- names(inits)
  # ageclass_nll(obsageclass = )
  age_res <- stats::optim(par = inits,fn = ageclass_nll,method = "L-BFGS-B",
                   lower = lower,
                   upper = upper,
                   max_a = max(s4t_ch$s4t_config$set_max_a),
                   mod_mat_a_beta = ageclass_data$mod_mat_a_beta,
                   obsageclass = ageclass_data$obsageclass,
                   ll = TRUE,
                   hessian = TRUE,
                   control = list(maxit = 500))

  # inits[1:16] <- c(-0.99,1.008,1.004,1.0005,-.82,0.099,0.098,.098,.099,.0986,.04,.03,-.14,
  # -.078,.02,-.11)

  # hes <- numDeriv::hessian(ageclass_nll,age_res$par,
  #              max_a = s4t_ch$ch_info$max_a,
  #              mod_mat_a_beta = ageclass_data$mod_mat_a_beta,
  #              obsageclass = ageclass_data$obsageclass,
  #              ll = TRUE)
  # solve(hes)

  estimated_parameters <- data.frame(parameter = names(age_res$par),
                                     estimate = age_res$par,
                                     std_error = sqrt(diag(solve(age_res$hessian))))
  estimated_parameters$z_value <- estimated_parameters$estimate / estimated_parameters$std_error
  estimated_parameters$lcl95 <- estimated_parameters$estimate - 1.96 * estimated_parameters$std_error
  estimated_parameters$ucl95 <- estimated_parameters$estimate + 1.96 * estimated_parameters$std_error


  s4t_ageclass_model <- list(estimated_parameters = estimated_parameters,
                             res = age_res,
                             call = match.call(),
                             AIC = age_res$value + 2 * length(age_res$par),
                             nll = age_res$value, k = length(age_res$par),
                             s4t_ch = s4t_ch)

  class(s4t_ageclass_model) <- "s4t_ageclass_model"

  return(s4t_ageclass_model)
}
