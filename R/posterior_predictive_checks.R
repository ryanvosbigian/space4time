


posterior_predict_knownage_s4t <- function(mod) {
  s4t_ch <- mod$fit$s4t_ch

  if (!is(mod,"s4t_cjs_rstan")) stop("Only implemented for objects of class s4t_cjs_rstan")



  init_rel <- s4t_ch$ch$m_matrix
  init_rel <- cbind(init_rel,init_rel= ifelse(init_rel[,"j"] == init_rel[,"r"],1,0))
  init_rel_only <- init_rel[init_rel[,"j"] == init_rel[,"r"],1:8]

  ka_init_rel_only <- init_rel_only[!is.na(init_rel_only[,"ageclass"]),]

  ka_init_rel_only


  recaptures <- matrix(data = NA, nrow = 0,ncol = ncol(ka_init_rel_only))
  colnames(recaptures) <- colnames(ka_init_rel_only)

  cohort_transitions <- mod$cohort_transitions
  mod$original_units$indices_theta_original

  for (i in 1:nrow(ka_init_rel_only)) {


    which_cohort <- which(cohort_transitions$a1 == ka_init_rel_only[i,"ageclass"] &
                         cohort_transitions$j == ka_init_rel_only[i,"j"] &
                         cohort_transitions$s == ka_init_rel_only[i,"s"] &
                         cohort_transitions$r == ka_init_rel_only[i,"r"] &
                         cohort_transitions$g == ka_init_rel_only[i,"g"])

    # rownames(cohort_transitions)[
    # )]
    #

  }







  init_rel_m_aux_probs <- apply(alk_pars,MARGIN = 1,FUN = function(pars) {
    alphas <- pars[1:ALPHA_col]
    deltas <- pars[(ALPHA_col+1):ALK_col]

    eta <- mod_mat_m_aux %*%  deltas

    ordered <- matrix(NA,nrow = nrow(eta),ncol = ALPHA_col)

    ordered[,1] <- 1 - stats::plogis(eta - alphas[1])

    for (i in 2:(ALPHA_col - 1)) {
      ordered[,i] <- stats::plogis(eta - alphas[i-1]) - stats::plogis(eta - alphas[i])
    }
    ordered[,ALPHA_col] <- stats::plogis(eta - alphas[i-1])

  })


}







#
# posterior_predict_s4t <- function(mod) {
#   s4t_ch <- mod$fit$s4t_ch
#
#   if (!is(mod,"s4t_cjs_rstan")) stop("Only implemented for objects of class s4t_cjs_rstan")
#
#   if (mod$fit$fixed_age == TRUE) stop("Not implemented for models fit with fixed_age == TRUE")
#
#
#   init_rel <- s4t_ch$ch$m_matrix
#   init_rel <- cbind(init_rel,init_rel= ifelse(init_rel[,"j"] == init_rel[,"r"],1,0))
#
#   pars_alk_alpha <- rownames(mod$estimated_parameters)[grepl("alk_par_alpha",rownames(mod$estimated_parameters))]
#
#   pars_alk_delta <- rownames(mod$estimated_parameters)[grepl("alk_par_delta",rownames(mod$estimated_parameters))]
#
#   mod$estimated_parameters
#
#   alk_par_alpha <- rstan::extract(mod$res,"alk_par_alpha",permuted = TRUE)[[1]]
#
#   alk_par_delta <- rstan::extract(mod$res,"alk_par_delta",permuted = TRUE)[[1]]
#
#   alk_par_alpha <- cbind(rep(-Inf,nrow(alk_par_alpha)),
#         alk_par_alpha,
#         rep(Inf,nrow(alk_par_alpha)))
#
#
#   alk_pars <- cbind(alk_par_alpha,alk_par_delta)
#   ALPHA_col <- ncol(alk_par_alpha)
#   ALK_col <- ncol(alk_pars)
#
#
#
#
#
#   init_rel_m_aux <- s4t_ch$ch$m_aux_df[init_rel[,"init_rel"]==1,]
#
#   info_ageclass <- ageclass_call(age_formula = mod$fit$ageclass_formula,
#                                  obs_aux = init_rel_m_aux,
#                                  ll = FALSE)
#
#   mod_mat_m_aux <-info_ageclass$mod_mat_a_beta[,-1]
#
#
#   init_rel_m_aux_probs <- apply(alk_pars,MARGIN = 1,FUN = function(pars) {
#     alphas <- pars[1:ALPHA_col]
#     deltas <- pars[(ALPHA_col+1):ALK_col]
#
#     eta <- mod_mat_m_aux %*%  deltas
#
#     ordered <- matrix(NA,nrow = nrow(eta),ncol = ALPHA_col)
#
#     ordered[,1] <- 1 - plogis(eta - alphas[1])
#
#     for (i in 2:(ALPHA_col - 1)) {
#       ordered[,i] <- plogis(eta - alphas[i-1]) - plogis(eta - alphas[i])
#     }
#     ordered[,ALPHA_col] <- plogis(eta - alphas[i-1])
#
#   })
#
#
# }
