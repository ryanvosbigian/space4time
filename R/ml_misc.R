
cohort_surv <- function(par,object) {

  # par <- object$res$par

  chobject <- object$s4t_ch

  n_batches <- chobject$ch_info$n_batches
  batches_list <- chobject$ch_info$batches_list
  max_a <- chobject$ch_info$max_a
  max_s_rel <- chobject$ch_info$max_s_rel
  max_t_recap <- chobject$ch_info$max_t_recap
  n_stations <- chobject$ch_info$n_stations
  recap_sites <- chobject$ch_info$recap_sites
  last_sites <- chobject$ch_info$last_sites
  recap_sites_not_last <- chobject$ch_info$recap_sites_not_last
  not_last_sites <- chobject$ch_info$not_last_sites
  sites_config <- chobject$user_defined$sites_config
  holdover_config <- chobject$user_defined$holdover_config
  mod_mat_theta <- object$fit$mod_mat_theta
  mod_mat_p <- object$fit$mod_mat_p

  indices_theta <- object$fit$indices_theta


  # mod_mat_theta <- object$fit$

  max_t <- max(c(max_s_rel,max_t_recap))



  Theta <- array(0,dim = c(a1 = max_a, a2 = max_a+1,
                           s = max_t,
                           j = n_stations,
                           k = n_stations,
                           b = n_batches))

  theta <- array(0,dim = c(a1 = max_a, a2 = max_a,
                           s = max_t,
                           j = n_stations,
                           k = n_stations,
                           b = n_batches))

  not_last_sites <- c(1:n_stations)[-last_sites]

  mod_mat_theta; indices_theta

  theta_params <- par[grepl("theta",names(par))]
  p_params <- par[grepl("p_",names(par))]

  for (i in 1:nrow(indices_theta)) {
    a1 <- indices_theta[i,"a1"]
    a2 <- indices_theta[i,"a2"]
    s <- indices_theta[i,"s"]
    j <- indices_theta[i,"j"]
    k <- indices_theta[i,"k"]
    b <- indices_theta[i,"b"]

    theta[a1,a2,s,j,k,b] <- stats::plogis(sum(mod_mat_theta[i,] * theta_params))

  }


  theta_inv <- 1-theta
  for (j in 1:n_stations) {
    k <- which(sites_config[j,] == 1)
    for (a1 in 1:max_a) {
      for (s in 1:max_s_rel[j]) {
        for (b in batches_list[[j]]) {

          Theta[a1,a1,s,j,k,b] <- theta[a1,a1,s,j,k,b]

          for (a2 in 2:max_a) {

            Theta[a1,a2,s,j,k,b] <- theta[a1,a2,s,j,k,b]*prod(theta_inv[a1,1:(a2-1),s,j,k,b])

          }
          Theta[a1,max_a+1,s,j,k,b] <- prod(theta_inv[a1,1:max_a,s,j,k,b])


        } # b

      } #  t
    } #  a1
    # }
  }




  cohort_surv <- vector()

  not_last_sites <- c(1:n_stations)[-last_sites]

  for (j in not_last_sites) {
    k <- which(sites_config[j,] == 1)
    for (b in batches_list[[j]]) {
      for (a1 in 1:max_a) {
        for (s in 1:max_s_rel[j]) { # should this be "s"?
          tmp_max_t <- min(c(max_t_recap[k],
                             s + (max_a - a1)))

          for (t in (s):(tmp_max_t)) {
            a2 <- a1 + t - s


            tmp_upperage <- min(c(t - s + a1, max_a)); tmp_upperage

            cohort_surv[paste0("cohort_surv[",a1,",",a2,",",s,",",t,",",j,",",k,",",b,"]")] <-
              stats::qlogis(Theta[a1,a2,s,j,k,b])

          } # t
        } # s
      } # a1
    } #  b
  } #j

  return(cohort_surv)
}

estimate_cohort_surv <- function(object) {

  ja <- numDeriv::jacobian(cohort_surv,x = object$res$par,object = object)
  vc <- solve(object$res$hessian)

  der_vc <- ja %*% vc %*% t(ja)

  ests <- cohort_surv(par = object$res$par,object=object)

  par_name = gsub("cohort_surv[[]|[]]","",names(ests))

  # paste0("cohort_surv[",a1,",",a2,",",s,",",t,",",j,",",k,",",b,"]")

  age1 <- stringr::str_split_i(par_name,",",1)
  age2 <- stringr::str_split_i(par_name,",",2)
  time_s <- stringr::str_split_i(par_name,",",3)
  time_t <- stringr::str_split_i(par_name,",",4)
  site_j <- stringr::str_split_i(par_name,",",5)
  site_k <- stringr::str_split_i(par_name,",",6)
  batch <- stringr::str_split_i(par_name,",",7)

  data.frame(parameters = names(ests),
             age1 = age1,
             age2 = age2,
             s_rel = time_s,
             t_rec = time_t,
             site_rel = site_j,
             site_rec = site_k,
             batch = batch,
             estimates_tr = round(stats::plogis(ests),3),
             lcl = round(stats::plogis(ests - 1.96*sqrt(diag(der_vc))),4),
             ucl = round(stats::plogis(ests + 1.96*sqrt(diag(der_vc))),4),
             estimates_logitscale = round(ests,3),
             se = round(sqrt(diag(der_vc)),3))

}


cdl_alg_marg <- function(par,
                         m_matrix,
                         l_matrix,
                         obs_aux,
                         max_t_recap,
                         max_s_rel,
                         sites,
                         site_path,
                         holdover_config,
                         set_min_a,
                         set_max_a,
                         n_batches,
                         batches_list,
                         mod_mat_theta,
                         indices_theta,
                         mod_mat_p,
                         indices_p_obs,
                         ageclass_data,
                         ageclassdat_L,
                         ageclassdat_M,
                         fixed_age,
                         fixed_ageclass_l,
                         fixed_ageclass_m) {

  if (!fixed_age) {
    # define alk params
    alk_par <- par[grepl("^a_",names(par))]

    if (alk_par["a_alpha_2"] < 0) return(NA) ## TEMPORARY

  }



  # needs to be at least zero
  # alk_par[2] <- exp(alk_par[2])

  # max_t
  max_t <- max(c(max_s_rel,max_t_recap))
  n_stations <- ncol(sites)
  recap_sites <- which(colSums(sites) > 0)
  last_sites <- unique(unlist(lapply(site_path,FUN = max)))

  recap_sites_not_last <- recap_sites[!(recap_sites %in% last_sites)]


  not_last_sites <- c(1:n_stations)[-last_sites]



  Theta <- array(0,dim = c(a1 = max(set_max_a), a2 = max(set_max_a)+1,
                           t = max_t,
                           j = n_stations,
                           k = n_stations,
                           b = n_batches))

  theta <- array(0,dim = c(a1 = max(set_max_a), a2 = max(set_max_a),
                           t = max_t,
                           j = n_stations,
                           k = n_stations,
                           b = n_batches))

  not_last_sites <- c(1:n_stations)[-last_sites]

  mod_mat_theta; indices_theta

  theta_params <- par[grepl("theta",names(par))]
  p_params <- par[grepl("p_",names(par))]

  for (i in 1:nrow(indices_theta)) {
    a1 <- indices_theta[i,"a1"]
    a2 <- indices_theta[i,"a2"]
    s <- indices_theta[i,"s"]
    j <- indices_theta[i,"j"]
    k <- indices_theta[i,"k"]
    b <- indices_theta[i,"b"]

    theta[a1,a2,s,j,k,b] <- stats::plogis(sum(mod_mat_theta[i,] * theta_params))

  }

  p_obs <- array(0,dim = c(a1 = max(set_max_a),
                           a2 = max(set_max_a),
                           t = max(indices_p_obs[,"t"]),
                           k = n_stations,
                           b = n_batches))


  for (i in 1:nrow(indices_p_obs)) {
    a1 <- indices_p_obs[i,"a1"]
    a2 <- indices_p_obs[i,"a2"]
    t <- indices_p_obs[i,"t"]
    k <- indices_p_obs[i,"k"]
    b <- indices_p_obs[i,"b"]

    p_obs[a1,a2,t,k,b] <- stats::plogis(sum(mod_mat_p[i,] * p_params))

  }

  p_obs[,,,last_sites,] <- 1



  theta_inv <- 1-theta
  ## reconsider these loops for j and k. Should determine k based on sites matrix
  for (j in 1:n_stations) {
    k <- which(sites[j,] == 1)
    # for (k in 2:(length(max_t_recap))) {
    for (a1 in 1:set_max_a[j]) {
      for (t in 1:max_s_rel[j]) {
        for (b in batches_list[[j]]) {


          Theta[a1,a1,t,j,k,b] <- theta[a1,a1,t,j,k,b]

          for (a2 in 2:max(set_max_a)) {

            Theta[a1,a2,t,j,k,b] <- theta[a1,a2,t,j,k,b]*prod(theta_inv[a1,1:(a2-1),t,j,k,b])


          }

          Theta[a1,max(set_max_a)+1,t,j,k,b] <- prod(theta_inv[a1,1:max(set_max_a),t,j,k,b])

        } # b

      } #  t
    } #  a1
    # }
  }



  lambda_array <- array(0,dim = c(j = n_stations,
                                  k = n_stations,
                                  s = max_t,
                                  t = max_t,
                                  b = n_batches,
                                  a1 = max(set_max_a)))

  not_last_sites <- c(1:n_stations)[-last_sites]

  for (j in not_last_sites) {
    k <- which(sites[j,]==1) #
    for (b in batches_list[[j]]) {
      for (s in 1:max_s_rel[j]) {
        for (a1 in 1:set_max_a[j]) {

          tmp_max_t <- min(c(max_t_recap[k],s + (set_max_a[k] - a1)))
          for (t in s:tmp_max_t) {
            a2 <- a1 + t - s
            lambda_array[j,k,s,t,b,a1] <- Theta[a1,a2,s,j,k,b]
          } # t
        } # a1
      } # s
    } # b
  } # j

  # Determine sites that have more than one site to be recaptured
  site_path_lengths <- unlist(lapply(site_path,length))
  site_path_length3 <- which(unlist(lapply(site_path,length)) >= 3)


  for (j in site_path_length3) {
    tmp_ks <- site_path[[j]]
    for (n.k in 3:site_path_lengths[j]) {
      k <- tmp_ks[n.k]
      k_min1 <- tmp_ks[n.k - 1]
      # for (k in tmp_ks[3:n_stations]) {
      for (b in batches_list[[j]]) {
        # use max_s_rel??
        for (s in 1:(max_s_rel[j])) {
          # use max_s_rel??
          for (a1 in 1:set_max_a[j]) { # a is age at time s
            tmp_max_t <- min(c(max_t_recap[k],s + (set_max_a[k] - a1)))
            for (t in (s):(tmp_max_t)) {
              a2 <- a1 + t - s

              # CAN JUST DO THE LOOP FOR AGE BEFORE THE LOOP FOR T
              if (t - s > set_max_a[k]- a1) stop()

              # tmp_upperage <- ifelse(t - s + a1 > max_a, max_a, t - s + a1); tmp_upperage

              tmp_upperage <- min(c(t - s + a1, set_max_a[k])); tmp_upperage

              # if (t - s + a1 > max_a) stop()

              lambda_array[j,k,s,t,b,a1] <- sum(lambda_array[j,k_min1,s,s:t,b,a1]*
                                                  (1-flex_diag(p_obs[a1,a1:tmp_upperage,s:t,k_min1,b]))*
                                                  flex_diag(lambda_array[k_min1,k,s:t,t,b,a1:tmp_upperage]))
              # if (sum(is.na(lambda_array[j,k,s,t,b,a1])) > 0) stop()

            } # t
          } # a1
        } # s
      } # b
    } # k
  }

  chi_array <- array(0,dim = c(j = n_stations,s = max_t,b = n_batches,a = max(set_max_a)))

  chi_array[n_stations,,,] <- 1


  not_last_sites_rev <- rev(c(1:n_stations)[-last_sites])

  for (j in not_last_sites_rev) {
    k <- which(sites[j,] == 1)

    for (b in batches_list[[j]]) {
      for (s in 1:max_s_rel[j]) {
        for (a1 in 1:set_max_a[j]) {

          tmp_max_t <- min(c(max_t_recap[k],
                             s + (set_max_a[k] - a1)))

          tmp_maxage <- min(c(set_max_a[k],tmp_max_t - s + a1))
          max(set_max_a)

          chi_array[j,s,b,a1] <- Theta[a1,max(set_max_a)+1,s,j,k,b] +
            sum(Theta[a1,a1:(tmp_maxage),s,j,k,b]  * (1 - flex_diag(p_obs[a1,a1:tmp_maxage,s:tmp_max_t,k,b])) *
                  flex_diag(chi_array[k,s:tmp_max_t,b,a1:(tmp_maxage)]))

        }
      }
    }

  }

  if (fixed_age) {
    pi_L <- fixed_ageclass_l
  } else {
    pi_L <- ageclass_nll(par = alk_par,
                         max_a = max(set_max_a),
                         mod_mat_a_beta = ageclassdat_L$mod_mat_a_beta,
                         ll = FALSE)

  }



  # known ages:
  knownage_l <- which(!is.na(ageclassdat_L$obsageclass))
  unknownage_l <- which(is.na(ageclassdat_L$obsageclass))

  ll <- 0
  ll_vec <- vector(length = nrow(l_matrix))

  for (i in knownage_l) {
    diff_age <- (l_matrix[i,"s"] - l_matrix[i,"obs_time"]); diff_age
    k_age <- diff_age + ageclassdat_L$obsageclass[i]

    lik1 <- chi_array[l_matrix[i,1],l_matrix[i,2],l_matrix[i,3],k_age]
    ll_vec[i] <- - log(lik1) * l_matrix[i,"n"]


  }




  for (i in unknownage_l) { # 1:nrow(l_matrix)
    # l_matrix[i,]


    # l_matrix %>% as.data.frame() %>% mutate(i = 1:n()) %>% filter(k == 2,t == 1)

    # adjust depending on difference in age
    diff_age <- (l_matrix[i,"s"] - l_matrix[i,"obs_time"]); diff_age


    pi_age2 <- c(rep(0,diff_age),pi_L[i,(1):(max(set_max_a) - diff_age)]); pi_age2
    pi_age3 <- pi_age2/sum(pi_age2) # prob wrong



    lik1 <- chi_array[l_matrix[i,1],l_matrix[i,2],l_matrix[i,3],]
    lik2 <- lik1 * pi_age3
    ll_vec[i] <- - log(sum(lik2)) * l_matrix[i,"n"]

  }

  if (any(is.nan(ll_vec))) {
    i <- which(is.nan(ll_vec) | is.infinite(ll_vec))
    print("l_matrix")
    print(l_matrix[i,])
    print(pi_L[i,])
    print(ll_vec[i])
    print(par)
    print(i)
    stop()
  }

  ll <- ll + sum(ll_vec)


  if (fixed_age) {
    pi_M <- fixed_ageclass_m
  } else {
    pi_M <- ageclass_nll(par = alk_par,
                         max_a = max(set_max_a),
                         mod_mat_a_beta = ageclassdat_M$mod_mat_a_beta,
                         ll = FALSE)

  }




  # known ages:
  knownage_m <- which(!is.na(ageclassdat_M$obsageclass))
  unknownage_m <- which(is.na(ageclassdat_M$obsageclass))

  ll_vec <- vector(length = nrow(m_matrix))

  for (i in knownage_m) {
    # m_matrix[i,]

    diff_age <- (m_matrix[i,"s"] - m_matrix[i,"obs_time"]); diff_age
    diff_age_t <- (m_matrix[i,"t"] - m_matrix[i,"obs_time"]); diff_age_t

    k <- m_matrix[i,"k"]; k


    ## CHANGING THIS TOO!!!
    tmp_p_obs1 <- flex_diag(p_obs[1:(set_max_a[k]-diff_age_t),(1+diff_age_t):set_max_a[k],m_matrix[i,"t"],m_matrix[i,"k"],m_matrix[i,"b"]])
    tmp_p_obs2 <- c(tmp_p_obs1,rep(0,diff_age_t),rep(0,max(set_max_a) - set_max_a[k]))


    j_age <- diff_age + ageclassdat_M$obsageclass[i]


    # tmp_p_obs2[k_age]

    lik1 <- tmp_p_obs2[j_age]*
      lambda_array[m_matrix[i,1],m_matrix[i,2],m_matrix[i,3],m_matrix[i,4],m_matrix[i,5],j_age]


    #
    ll_vec[i] <- - log(lik1)*m_matrix[i,"n"]
  }


  for (i in unknownage_m) { # 1:nrow(m_matrix)
    # m_matrix[i,]

    # adjust depending on difference in age
    diff_age <- (m_matrix[i,"s"] - m_matrix[i,"obs_time"]); diff_age

    k <- m_matrix[i,"k"]

    # pi_age1 <- pi_M[i,]
    # pi_age2 <- c(rep(0,diff_age),pi_M[i,(diff_age+1):max(set_max_a)])

    ## CHANGING THIS!!!?!!!
    pi_age2 <- c(rep(0,diff_age),pi_M[i,(1):(max(set_max_a) - diff_age)])
    pi_age3 <- pi_age2/sum(pi_age2) # prob wrong

    diff_age_t <- (m_matrix[i,"t"] - m_matrix[i,"obs_time"]); diff_age_t


    ## CHANGING THIS TOO!!!
    tmp_p_obs1 <- flex_diag(p_obs[1:(set_max_a[k]-diff_age_t),(1+diff_age_t):set_max_a[k],m_matrix[i,"t"],m_matrix[i,"k"],m_matrix[i,"b"]])
    tmp_p_obs2 <- c(tmp_p_obs1,rep(0,diff_age_t),rep(0,max(set_max_a) - set_max_a[k]))


    # tmp_p_obs1 <- flex_diag(p_obs[1:(set_max_a[k]-diff_age_t),(1+diff_age_t):set_max_a[k],m_matrix[i,"t"],m_matrix[i,"k"],m_matrix[i,"b"]])
    # tmp_p_obs2 <- c(rep(0,diff_age),tmp_p_obs1,rep(0,diff_age_t-diff_age))



    lik1 <- tmp_p_obs2 *
      lambda_array[m_matrix[i,1],m_matrix[i,2],m_matrix[i,3],m_matrix[i,4],m_matrix[i,5],]



    lik2 <- lik1 * pi_age3

    ll_vec[i] <- - log(sum(lik2))*m_matrix[i,"n"]

    # some different debugging stuff

    # if (is.nan(log(sum(lik2))) | is.infinite(log(sum(lik2)))) {
    #   print("m_matrix")
    #   print(m_matrix[i,])
    #   print(lik2)
    #   print(pi_age)
    #   print(par)
    #   print(i)numDeriv
    #   stop()
    # }

  }

  if (any(is.nan(ll_vec)) | any(is.infinite(ll_vec))) {
    i <- which(is.nan(ll_vec) | is.infinite(ll_vec))
    print("m_matrix")
    print(m_matrix[i,])
    print(ll_vec[i])
    print(pi_M[i,])
    print(par)
    print(i)
    stop()
  }

  ll <- ll + sum(ll_vec)

  if (fixed_age) {
    alk_ll <- 0
  } else {
    alk_ll <- ageclass_nll(par = alk_par,
                           max_a = max(set_max_a),
                           mod_mat_a_beta = ageclass_data$mod_mat_a_beta,
                           obsageclass = ageclass_data$obsageclass,
                           ll = TRUE)
  }




  joint_nll <- ll + alk_ll

  if (is.finite(joint_nll) == FALSE) {
    print("NLL")
    print(joint_nll)
    # print(A_array)
    print(par)
  }

  return(joint_nll)


}
