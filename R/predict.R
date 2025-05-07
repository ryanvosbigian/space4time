



theta_to_Theta <- function(theta) {

  max_a <- dim(theta)[1]
  n_stations <- dim(theta)[4]
  max_s <- dim(theta)[3]
  n_batches <- dim(theta)[5]
  n_groups <- dim(theta)[6]

  theta_inv <- 1-theta
  ## reconsider these loops for j and k. Should determine k based on sites matrix
  for (g in 1:n_groups) {
    for (j in 1:n_stations) {
      # k <- which(sites_config[j,] == 1)
      # for (k in 2:(length(max_t_recap))) {
      for (a1 in 1:max_a) {
        for (s in 1:max_s) {
          for (b in 1:n_batches) {


            Theta[a1,a1,s,j,b,g] <- theta[a1,a1,s,j,b,g]

            for (a2 in 2:max_a) {

              Theta[a1,a2,s,j,b,g] <- theta[a1,a2,s,j,b,g]*prod(theta_inv[a1,1:(a2-1),s,j,b,g])


            }

            Theta[a1,max_a+1,s,j,b,g] <- prod(theta_inv[a1,1:max_a,s,j,b,g])

          } # b

        } #  t
      } #  a1
      # }
    }
  } # g



}



theta_to_cohortSurv <- function(theta,indices_cohort_surv) {

  Theta <- array(0,dim = dim(theta) + c(0,1,0,0,0,0))

  max_a <- dim(theta)[1]
  n_stations <- dim(theta)[4]
  max_s <- dim(theta)[3]
  n_batches <- dim(theta)[5]
  n_groups <- dim(theta)[6]

  theta_inv <- 1-theta
  ## reconsider these loops for j and k. Should determine k based on sites matrix
  for (g in 1:n_groups) {
    for (j in 1:n_stations) {
      # k <- which(sites_config[j,] == 1)
      # for (k in 2:(length(max_t_recap))) {
      for (a1 in 1:max_a) {
        for (s in 1:max_s) {
          for (b in 1:n_batches) {


            Theta[a1,a1,s,j,b,g] <- theta[a1,a1,s,j,b,g]

            for (a2 in 2:max(max_a)) {

              Theta[a1,a2,s,j,b,g] <- theta[a1,a2,s,j,b,g]*prod(theta_inv[a1,1:(a2-1),s,j,b,g])


            }

            Theta[a1,max_a+1,s,j,b,g] <- prod(theta_inv[a1,1:max_a,s,j,b,g])

          } # b

        } #  t
      } #  a1
      # }
    }
  } # g


  cohort_surv <- vector(length = nrow(indices_cohort_surv))

  for (i in 1:nrow(indices_cohort_surv)) {
    cohort_surv[i] <- Theta[indices_cohort_surv[i,1],
                            indices_cohort_surv[i,2],
                            indices_cohort_surv[i,3],
                            indices_cohort_surv[i,5],
                            indices_cohort_surv[i,7],
                            indices_cohort_surv[i,8]]
  }

  names(cohort_surv) <- apply(indices_cohort_surv,MARGIN = 1,FUN = function(x) paste0("cohort[",(paste0(x,collapse = ",")),"]"))
  return(cohort_surv)
}

populate_theta <- function(theta,indices_theta,theta_vals) {
  for (i in 1:nrow(indices_theta)) {
    theta[indices_theta[i,1],
          indices_theta[i,2],
          indices_theta[i,3],
          indices_theta[i,5],
          indices_theta[i,7],
          indices_theta[i,8]] <- theta_vals[i]
  }
  return(theta)
}


predictTheta.s4t_cjs_rstan <- function(object,newdata = NULL,lcl = 0.025,
                                       ucl = 0.975,...) {
  (object$res@sim$samples)

  theta_pars <- object$original_units$compare_parnames[grepl("theta_",object$original_units$compare_parnames[,1]),1]

  samps <- rstan::extract(object$res,pars = theta_pars)
  # samps <- rstan::as.array(object$res)
  theta_samps <- cbind(sapply(samps,FUN = function(x) x))
  # class(theta_mat)

  if (is.null(newdata)) newdata <- object$fit$input_data$mod_mat_theta

  theta_samps <- apply(theta_samps,MARGIN = 1,function(x) stats::plogis(newdata %*% x))
  # populate_theta

  indices_cohort_surv <- object$fit$input_data$indices_cohort_surv
  indices_theta <- object$fit$input_data$indices_theta

  theta_dims <- apply(indices_theta,2,max)[c("a1","a2","s","j","b","g")]

  indices_theta <- object$fit$input_data$indices_theta

  theta <- array(data = 0, dim = theta_dims)

  cohort_surv_samps <- apply(theta_samps,MARGIN = 2,FUN = function(x) theta_to_cohortSurv(populate_theta(theta= theta,
                                                                indices_theta = indices_theta,
                                                                theta_vals = x),
                                                                indices_cohort_surv = indices_cohort_surv
                                                                  ))


  cohort_surv_samps

  cohort_surv_mean <- apply(cohort_surv_samps,MARGIN = 1,mean)
  cohort_surv_lcl <- apply(cohort_surv_samps,MARGIN = 1,function(x) stats::quantile(x,prob = lcl))
  cohort_surv_ucl <- apply(cohort_surv_samps,MARGIN = 1,function(x) stats::quantile(x,prob = ucl))

  # theta_mean <- apply(theta_samps,MARGIN = 2,mean)
  # theta_lcl <- apply(theta_samps,MARGIN = 2,function(x) quantile(x,prob = lcl))
  # theta_ucl <- apply(theta_samps,MARGIN = 2,function(x) quantile(x,prob = ucl))

  # indices_theta_original <- object$original_units$indices_theta_original

  preds <- cbind(indices_cohort_surv,
                 cohort_surv_mean = cohort_surv_mean,
                 cohort_surv_lcl = cohort_surv_lcl,
                 cohort_surv_ucl = cohort_surv_ucl)

  return(preds)

}


#' Predict method for fit_ageclass fits
#'
#' add description
#'
#' @export
#'
#' @param object a `s4t_ageclass_model` object
#' @param newdata add description.
#' @param type a `character` of one of the three option `prob`, `class`,
#'    or `cum.prob`
#' @param ... not used
#' @returns depends on the type.
predict.s4t_ageclass_model <- function(object, newdata,
                                 type = c("prob","class","cum.prob"),
                                 ...) {


  type <- match.arg(type)

  stopifnot(is(object,"s4t_ageclass_model"))

  if (missing(newdata)) {
    newdata <- object$s4t_ch$ch$obs_aux
  }

  info_ageclass <- ageclass_call(age_formula = object$call$age_formula,
                                 obs_aux = newdata,
                                 ll = FALSE)

  prob <- ageclass_nll(par = object$res$par,
                       max_a = max(object$s4t_ch$ch_info$set_max_a),
                       mod_mat_a_beta = info_ageclass$mod_mat_a_beta,
                       ll = FALSE)

  min_a <- object$s4t_ch$ch_info$observed_relative_min_max$min_obs_age
  max_a <-  object$s4t_ch$ch_info$observed_relative_min_max$max_obs_age

  age_cols <- paste0("Age",min_a:(ncol(prob) + min_a - 1))

  colnames(prob) = age_cols

  if (type == "prob") {
    return(prob)
  } else if (type == "class") {
    pred_age <- apply(prob,MARGIN = 1,which.max) + min_a - 1
    return(pred_age)
  } else {
    # type == cum.prob
    cum.prob <- t(apply(prob,1,function(x) cumsum(x)))
    return(cum.prob)
  }





}
