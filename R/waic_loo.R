


#' Extract pointwise log-likelihood from a fitted s4t_cjs_rstan model
#'
#' Function for extracting and de-marginalizing the pointwise log-likelihood matrix
#'     from a s4t_cjs_rstan object. See loo::extract_log_lik.
#'
#' @export
#'
#' @importFrom loo extract_log_lik
#'
#' @param mod a `s4t_cjs_rstan` object
#' @returns an `S` by `N` matrix of post-warmup extracted draws, where `S` is
#'     the size of the posterior sample, and `N` is the number of unique observations.
#'
extract_log_lik_s4t <- function(mod) {
  if (!is(mod,"s4t_cjs_rstan")) stop("not an s4t_cjs_rstan object")

  if (!is.null(mod$fit$fixed_age)){
    fixed_age = TRUE
  } else{
    fixed_age = FALSE
  }

  loglik_m = loo::extract_log_lik(mod$res)

  if (fixed_age) {
    n_a <- 0
  } else {
    n_a <- mod$fit$input_data$N_obsageclass
  }

  n_l <- sum(mod$fit$input_data$l_matrix[,6])
  n_m <- sum(mod$fit$input_data$m_matrix[,8])

  n_l_marg <- nrow(mod$fit$input_data$l_matrix)
  n_m_marg <- nrow(mod$fit$input_data$m_matrix)

  new_loglik_m <- matrix(NA,nrow = nrow(loglik_m),n_a + n_l + n_m)

  if (fixed_age == FALSE) {
    new_loglik_m[,1:n_a] <- loglik_m[,1:n_a]
  }



  tmp_l_mat_ll <- loglik_m[,(n_a+1):(n_a + n_l_marg)]
  tmp_m_mat_ll <- loglik_m[,(n_a + n_l_marg + 1):(n_a + n_l_marg + n_m_marg)]

  n_obs_l_mat <- mod$fit$input_data$l_matrix[,6]
  tmp_unmarg_l_mat_ll <- do.call(cbind,
                                 lapply(1:nrow(mod$fit$input_data$l_matrix),
                                        FUN = function(i) replicate(n_obs_l_mat[i],
                                                                    tmp_l_mat_ll[,i]/n_obs_l_mat[i])

                                 )
  )

  # dim(tmp_unmarg_l_mat_ll)

  n_obs_m_mat <- mod$fit$input_data$m_matrix[,8]
  tmp_unmarg_m_mat_ll <- do.call(cbind,
                                 lapply(1:nrow(mod$fit$input_data$m_matrix),
                                        FUN = function(i) replicate(n_obs_m_mat[i],
                                                                    tmp_m_mat_ll[,i]/n_obs_m_mat[i])

                                 )
  )
  # dim(tmp_unmarg_m_mat_ll)


  new_loglik_m[,(n_a+1):(n_a+ n_l)] <- tmp_unmarg_l_mat_ll
  new_loglik_m[,(n_a+ n_l + 1):(n_a+ n_l + n_m)]  <- tmp_unmarg_m_mat_ll


  return(new_loglik_m)

}






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
#' @returns a named list with class `c("psis_loo", "loo")`. See `?loo::loo`
#'
loo.s4t_cjs_rstan <- function(x, ...) {

  log_lik_1 <- extract_log_lik_s4t(x)
  # loo::waic(log_lik_1)
  loo1 <- loo::loo(log_lik_1, ...)
  return(loo1)
}

#' Widely applicable information criterion (WAIC)
#'
#' The waic() methods can be used to compute WAIC from the pointwise log-likelihood. See `loo::waic()`
#' @export
#'
#' @importFrom loo waic
#' @importFrom loo extract_log_lik
#'
#' @param x a `s4t_cjs_rstan` object
#' @param ... not used (for generic consistency)
#' @returns A named list (of class `c("waic", "loo")`). See `?loo::waic`.
#'
waic.s4t_cjs_rstan <- function(x, ...) {

  log_lik_1 <- extract_log_lik_s4t(x)
  loo::waic(log_lik_1)
}
