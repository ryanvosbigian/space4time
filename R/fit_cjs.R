

#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom RcppParallel CxxFlags
#' @import methods
#' @import Rcpp



format_s4t_cjs <- function(p_formula,
                           theta_formula,
                           ageclass_formula,
                           cov_p = NULL, cov_theta = NULL,groups = NULL,
                           s4t_ch) {
  holdover_config <- s4t_ch$user_defined$holdover_config
  sites <- sites_config <- s4t_ch$user_defined$sites_config

  recap_sites <- s4t_ch$ch_info$recap_sites #which(colSums(sites) > 0)
  last_sites <- s4t_ch$ch_info$last_sites # unique(unlist(lapply(site_path,FUN = max)))

  n_stations <- s4t_ch$ch_info$n_stations

  n_batches <- s4t_ch$ch_info$n_batches
  batches_list <- s4t_ch$ch_info$batches_list

  max_s_rel <- as.integer(s4t_ch$ch_info$max_s_rel)
  max_t_recap <- as.integer(s4t_ch$ch_info$max_t_recap)
  set_max_a <- s4t_ch$ch_info$set_max_a
  set_min_a <- s4t_ch$ch_info$set_min_a
  first_obs <- s4t_ch$ch_info$first_obs
  first_sites <- s4t_ch$ch_info$first_sites
  site_path <- s4t_ch$ch_info$site_path

  min_ageclass_mat <- s4t_ch$ch_info$min_ageclass_mat
  max_ageclass_mat <- s4t_ch$ch_info$max_ageclass_mat

  ###
  # groups
  if (!is.null(groups)) {
    if (length(setdiff(groups,colnames(s4t_ch$ch$obs_aux))) > 0) {
      stop(paste0("groups were not found in obs_aux: ",paste0(setdiff(groups,colnames(s4t_ch$ch$obs_aux)),collapse = ", ")))
    }

    # the below needs to go after the covariates are added in.
    if (length(intersect(groups,"")) > 0) stop("group names are also in covariates")

    df_groups <- as.data.frame(as.data.frame(s4t_ch$ch$obs_aux)[,groups])

    if (length(groups) == 1) colnames(df_groups) <- groups

    distinct_groups_df <- dplyr::distinct_all(df_groups)

    distinct_groups_df$g <- as.factor(1:nrow(distinct_groups_df))

    N_groups <- nrow(distinct_groups_df)

    if (N_groups > 20) message(paste0("Number of distinct groups may be too large, N_groups = ",N_groups))

    # adding group to obs_aux
    suppressMessages(df_groups <- dplyr::left_join(df_groups,distinct_groups_df))

    group_id_obs_aux <- df_groups$g

    # adding group to l_matrix
    suppressMessages(tmp_groups <- dplyr::left_join((s4t_ch$ch$l_aux_df),distinct_groups_df))

    s4t_ch$ch$l_matrix[,"g"] <- tmp_groups$g

    # adding group to m_matrix
    suppressMessages(tmp_groups <- dplyr::left_join(as.data.frame(s4t_ch$ch$m_aux_df),distinct_groups_df))

    s4t_ch$ch$m_matrix[,"g"] <- tmp_groups$g

  } else {
    N_groups <- 1
    group_id <- as.factor(rep(1,nrow(s4t_ch$ch$obs_aux)))

  } # end if statement for groups


  ###


  recap_sites_not_last <- s4t_ch$ch_info$recap_sites_not_last


  # 7 columns because first two are a1 and a2. Third is t. Fourth is j, fifth is k, sixth is batch, seventh is group
  indices_theta <- matrix(NA, nrow = 0,ncol = 8)
  colnames(indices_theta) <- c("a1","a2","s","t","j","k","b","g")


  for (g in 1:N_groups) {
    for (j in 1:(n_stations-1)) {
      for (k in which(sites[j,]==1)) { # doesn't need to be a loop

        # sites[,k]
        # making an assumption that transitions between sites are the same
        # for all batches.
        for (b in batches_list[[j]]) {
          if (holdover_config[j,k] == 0) {
            # min_a <- max(c(1,))
            # max_a +
            # max_s_rel[j]
            # s in 1:max_s_rel[j]
            for (s in 1:max_s_rel[j]) {
              t <- s
              # tmp_min_a <- set_min_a[j] # max(1,t - max_s_rel[b] + 1)
              # tmp_max_a <- set_max_a[j]

              tmp_min_a <- min_ageclass_mat[j,s]
              tmp_max_a <- max_ageclass_mat[j,s]

              if (is.na(tmp_min_a) | is.na(tmp_max_a)) {
                next()
                # if it is NA, then skip, because it is likely a
                # site-time combination at a release site when no individuals were released
              }

              # print(paste0(t," ",min_a))

              tmp_indices <- data.frame(a1 = tmp_min_a:tmp_max_a,
                                        a2 = tmp_min_a:tmp_max_a,
                                        s = s,
                                        t = t,
                                        j = j,
                                        k = k,
                                        b = b,
                                        g = g)

              indices_theta <- rbind(indices_theta, tmp_indices)
            }


          } else { # holdover_config[j,k] == 1
            # tmp_min_a <- set_min_a[j]

            # I think that starting from 1 is ok because of how i set the first
            # occasion/time to 1
            for (s in 1:max_s_rel[j]) {


              # tmp_max_a <- set_max_a[k]

              tmp_min_a <- min_ageclass_mat[j,s]
              tmp_max_a <- max_ageclass_mat[k,s]

              if (is.na(tmp_min_a) | is.na(tmp_max_a)) {
                next()
                # if it is NA, then skip, because it is likely a
                # site-time combination at a release site when no individuals were released
              }

              tmp_indices <- data.frame(a1 = rep_crosses(tmp_min_a,tmp_max_a),
                                        a2 = misc_ind_2(tmp_min_a,tmp_max_a),
                                        s = s,
                                        t = s + misc_ind_2(tmp_min_a,tmp_max_a) -
                                          rep_crosses(tmp_min_a,tmp_max_a),
                                        j = j,
                                        k = k,
                                        b = b,
                                        g = g)

              # drop transitions that exceed max_t_recap[k]
              tmp_recap_t <- tmp_indices$a2 - tmp_indices$a1 + tmp_indices$s
              tmp_keeprows <- tmp_recap_t <= max_t_recap[k]
              tmp_indices <- tmp_indices[tmp_keeprows,]


              indices_theta <- rbind(indices_theta, tmp_indices)
            }
          }
        } # b

      } # k
    } # j
  }



  indices_p_obs <- matrix(NA, nrow = 0,ncol = 8)
  colnames(indices_p_obs) <- c("a1","a2","s","t","j","k","b","g")

  # recap_sites[!(recap_sites %in% last_sites)]
  for (g in 1:N_groups) {
    for (k in recap_sites_not_last) {

      for (b in batches_list[[k]]) {
        max_t_recap[k]

        # Identify the site path for this release group (batch). Then,
        # obtain the previous site (j)
        j <- site_path[[b]][which(site_path[[b]] == k) - 1]

        # skip if this isn't a recap for this particular batch/release group
        if (length(j) == 0) next()

        for (s in 1:max_s_rel[j]) {
          if (sum(holdover_config[,k]) == 0) {

            # tmp_max_a <- set_max_a[k]
            # tmp_min_a <- set_min_a[j]

            tmp_min_a <- min_ageclass_mat[j,s]
            tmp_max_a <- max_ageclass_mat[k,s]

            if (is.na(tmp_min_a) | is.na(tmp_max_a)) {
              next()
              # if it is NA, then skip, because it is likely a
              # site-time combination at a release site when no individuals were released
            }

            tmp_indices <- data.frame(a1 = tmp_min_a:tmp_max_a,
                                      a2 = tmp_min_a:tmp_max_a,
                                      s = s,
                                      t = s,
                                      j = j,
                                      k = k,
                                      b=b,
                                      g=g)
            indices_p_obs <- rbind(indices_p_obs, tmp_indices)
          } else {
            # tmp_max_a <- set_max_a[k]
            # tmp_min_a <- set_min_a[j]

            tmp_min_a <- min_ageclass_mat[j,s]
            tmp_max_a <- max_ageclass_mat[k,s]

            if (is.na(tmp_min_a) | is.na(tmp_max_a)) {
              next()
              # if it is NA, then skip, because it is likely a
              # site-time combination at a release site when no individuals were released
            }

            tmp_indices <- suppressWarnings(data.frame(a1 = rep_crosses(tmp_min_a,tmp_max_a),
                                                       a2 = misc_ind_2(tmp_min_a,tmp_max_a),
                                                       s = s,
                                                       t = s + (misc_ind_2(tmp_min_a,tmp_max_a) -
                                                                  rep_crosses(tmp_min_a,tmp_max_a)),
                                                       j = j,
                                                       k = k,
                                                       b = b,
                                                       g = g))

            # drop rows that exceed the max recap time
            tmp_indices <- subset(tmp_indices, t <= max_t_recap[k])
            indices_p_obs <- rbind(indices_p_obs, tmp_indices)

          }

        }

      }

    }
  }





  # mf2 <- model.frame(~ x * z,data = data.frame(x = 1:10,y = 5:14,z = as.character(rep(c(1:3),length=10))))
  # mt2 <- attr(mf2,"terms")
  # x2 <- model.matrix(mt2, mf2, contrasts)
  tmp_indices_p_obs <- as.data.frame(apply(indices_p_obs,MARGIN = 2,as.factor))
  tmp_indices_theta <- as.data.frame(apply(indices_theta,MARGIN = 2,as.factor))


  # cov_p <- data.frame(t = as.character(1:4),cov = rnorm(4))

  if (!is.null(cov_p)) {
    if (is(cov_p,"list")) {

      for (i in 1:length(cov_p)) {

        match_col <- intersect(colnames(cov_p[[i]]),colnames(tmp_indices_p_obs))

        for (j in match_col) {
          cov_p[[i]][,j] <- as.factor(cov_p[[i]][,j])
        }



        tmp_indices_p_obs <- dplyr::left_join(tmp_indices_p_obs,cov_p[[i]])
      }


    }
    if (is(cov_p,"data.frame")) {
      match_col <- intersect(colnames(cov_p),colnames(tmp_indices_p_obs))

      for (j in match_col) {
        cov_p[,j] <- as.factor(cov_p[,j])
      }

      # cov_p[,match_col] <- as.data.frame(apply(cov_p[,match_col],MARGIN = 2,as.factor))# as.factor(cov_p[,match_col])

      tmp_indices_p_obs <- dplyr::left_join(tmp_indices_p_obs,cov_p)

    }

  }


  if (!is.null(groups)) {
    tmp_indices_p_obs <- dplyr::left_join(tmp_indices_p_obs,distinct_groups_df, by ="g")

  }


  ## now with thetas

  # cov_theta <- data.frame(t = as.character(1:4),cov = rnorm(4))


  if (!is.null(cov_theta)) {
    if (is(cov_theta,"list")) {

      for (i in 1:length(cov_theta)) {

        match_col <- intersect(colnames(cov_theta[[i]]),colnames(tmp_indices_theta))

        # cov_theta[[i]][,match_col] <- as.factor(cov_theta[[i]][,match_col])

        for (j in match_col) {
          cov_theta[[i]][,j] <- as.factor(cov_theta[[i]][,j])
        }

        tmp_indices_theta <- dplyr::left_join(tmp_indices_theta,cov_theta[[i]])
      }


    }
    if (is(cov_theta,"data.frame")) {
      match_col <- intersect(colnames(cov_theta),colnames(tmp_indices_theta))

      # cov_theta[,match_col] <- as.factor(cov_theta[,match_col])

      for (j in match_col) {
        cov_theta[,j] <- as.factor(cov_theta[,j])
      }


      tmp_indices_theta <- dplyr::left_join(tmp_indices_theta,cov_theta)

    }

  }

  if (!is.null(groups)) {
    tmp_indices_theta <- dplyr::left_join(tmp_indices_theta,distinct_groups_df, by ="g")
  }





  mf_p <- stats::model.frame(p_formula,data = tmp_indices_p_obs,
                      drop.unused.levels = TRUE)
  mt_p <- attr(mf_p,"terms")
  mod_mat_p <- stats::model.matrix(mt_p, mf_p)
  p_par <- rep(0,length = length(colnames(mod_mat_p))); names(p_par) <- paste0("p_",colnames(mod_mat_p))


  ## dropping columns that are all zero


  if (any(colSums(mod_mat_p) == 0)) {
    keep_cols <- which(colSums(mod_mat_p) != 0)

    mod_mat_p <- mod_mat_p[,keep_cols]
    p_par <- p_par[keep_cols]

  }

  ## dropping duplicate columns
  if (ncol(mod_mat_p) > 1) {
    combos <- utils::combn(x = 1:ncol(mod_mat_p),m = 2)
    aliasedcols <- apply(combos,MARGIN = 2,FUN = function(x) all(mod_mat_p[,x[1]] == mod_mat_p[,x[2]])
    )

    if (sum(aliasedcols) > 0) {
      drop_cols <- combos[2,aliasedcols]

      keep_cols <- c(1:ncol(mod_mat_p))[-drop_cols]

      mod_mat_p <- mod_mat_p[,keep_cols]
      p_par <- p_par[keep_cols]
    }
  }

  # remove linear dependent columns to make a full rank design matrix
  # following method from WeightIt package (make.full.rank fn)

  keep <- rep.int(TRUE, ncol(mod_mat_p))
  mat_qr <- qr(mod_mat_p)
  keep[mat_qr$pivot[-seq(mat_qr$rank)]] <- FALSE

  mod_mat_p <- mod_mat_p[,keep, drop = FALSE]
  p_par <- p_par[keep]

  if (qr(mod_mat_p)$rank != ncol(mod_mat_p)) message("Model matrix for p is rank deficient")


  # x2[1,] * p_par

  mf_theta <- stats::model.frame(theta_formula,data = tmp_indices_theta,
                          drop.unused.levels = TRUE)
  mt_theta <- attr(mf_theta,"terms")
  mod_mat_theta <- stats::model.matrix(mt_theta, mf_theta)
  theta_par <- rep(0,length = length(colnames(mod_mat_theta))); names(theta_par) <- paste0("theta_",colnames(mod_mat_theta))



  ## dropping columns that are all zero


  if (any(colSums(mod_mat_theta) == 0)) {
    keep_cols <- which(colSums(mod_mat_theta) != 0)

    mod_mat_theta <- mod_mat_theta[,keep_cols]
    theta_par <- theta_par[keep_cols]

  }

  ## dropping aliased/duplicate columns

  if (ncol(mod_mat_theta) > 1) {
    combos <- utils::combn(x = 1:ncol(mod_mat_theta),m = 2)
    aliasedcols <- apply(combos,MARGIN = 2,FUN = function(x) all(mod_mat_theta[,x[1]] == mod_mat_theta[,x[2]])
    )

    if (sum(aliasedcols) > 0) {
      drop_cols <- combos[2,aliasedcols]

      keep_cols <- c(1:ncol(mod_mat_theta))[-drop_cols]

      mod_mat_theta <- mod_mat_theta[,keep_cols]
      theta_par <- theta_par[keep_cols]
    }
  }


  # remove linear dependent columns to make a full rank design matrix
  # following method from WeightIt package (make.full.rank fn)

  keep <- rep.int(TRUE, ncol(mod_mat_theta))
  mat_qr <- qr(mod_mat_theta)
  keep[mat_qr$pivot[-seq(mat_qr$rank)]] <- FALSE

  mod_mat_theta <- mod_mat_theta[,keep, drop = FALSE]
  theta_par <- theta_par[keep]


  if (qr(mod_mat_theta)$rank != ncol(mod_mat_theta)) message("Model matrix for theta is rank deficient")



  ### age length stuff:

  ageclassdat_L <- ageclass_call(age_formula=ageclass_formula,
                                 obs_aux = s4t_ch$ch$l_aux_df
                                 # max_a = s4t_ch$ch_info$max_a
  )

  ageclassdat_M <- ageclass_call(age_formula=ageclass_formula,
                                 obs_aux = s4t_ch$ch$m_aux_df
                                 # max_a = s4t_ch$ch_info$max_a
  )


  #
  #
  # not marginalizing this portion
  ageclass_data <- ageclass_call(age_formula=ageclass_formula,
                                 obs_aux = s4t_ch$ch$obs_aux[stats::complete.cases(s4t_ch$ch$obs_aux),]
                                 #na.omit(s4t_ch$ch$obs_aux),
                                 # max_a = s4t_ch$ch_info$max_a
  )
  inits <- ageclass_data$age_par



  ## ALSO NEED ALK PARAMS
  par = inits <- c(ageclass_data$age_par,p_par,theta_par)

  lower <- ifelse(grepl("a_alpha_",names(inits)) &
                    !grepl("^a_alpha_1$",names(inits)),
                  0,-Inf)

  upper <- rep(Inf,length(inits))

  names(lower) <- names(upper) <- names(inits)

  mod_mat_theta; indices_theta

  return(list(p_par = p_par,
              theta_par = theta_par,
              age_par = ageclass_data$age_par,
              m_matrix = s4t_ch$ch$m_matrix,
              l_matrix = s4t_ch$ch$l_matrix,
              m_aux_df = s4t_ch$ch$m_aux_df,
              l_aux_df = s4t_ch$ch$l_aux_df,
              obs_aux = s4t_ch$ch$obs_aux,
              max_t_recap = s4t_ch$ch_info$max_t_recap,
              max_s_rel = s4t_ch$ch_info$max_s_rel,
              sites = s4t_ch$user_defined$sites_config,
              site_path=s4t_ch$ch_info$site_path,
              holdover_config = s4t_ch$user_defined$holdover_config,
              set_min_a = s4t_ch$ch_info$set_min_a,
              set_max_a = s4t_ch$ch_info$set_max_a,
              n_batches  = s4t_ch$ch_info$n_batches,
              N_groups=N_groups,
              batches_list = s4t_ch$ch_info$batches_list,
              mod_mat_theta = mod_mat_theta,
              indices_theta = indices_theta,
              mod_mat_p = mod_mat_p,
              indices_p_obs = indices_p_obs,
              ageclass_data = ageclass_data,
              ageclassdat_L = ageclassdat_L,
              ageclassdat_M = ageclassdat_M
              # min_ageclass_mat = min_ageclass_mat,
              # max_ageclass_mat = max_ageclass_mat
  ))
}




fit_s4t_cjs_ml <- function(p_formula,theta_formula,
                             ageclass_formula,
                             cov_p = NULL, cov_theta = NULL,
                             groups = NULL,
                             s4t_ch,
                             ndeps = 1e-3,
                             lmm = 5,
                             maxit = 500,
                             fixed_age = FALSE) {

  format_cjs <- format_s4t_cjs(p_formula =p_formula, theta_formula = theta_formula,
                               ageclass_formula = ageclass_formula,
                               cov_p = cov_p,cov_theta = cov_theta,groups = groups,s4t_ch = s4t_ch)


  max_t_recap <- format_cjs$max_t_recap
  set_min_a <- format_cjs$set_min_a
  set_max_a <- format_cjs$set_max_a
  tmp_obs_aux <-format_cjs$obs_aux[!is.na(format_cjs$obs_aux[,"ageclass"]),]
  ageclassdat <- ageclass_call(age_formula = ageclass_formula,
                               obs_aux = tmp_obs_aux
                               # max_a = format_cjs$max_a,ll = TRUE
  )
  ageclassdat_M <- ageclass_call(age_formula = ageclass_formula,
                                 obs_aux = format_cjs$m_aux_df,
                                 # max_a = format_cjs$max_a,
                                 ll = TRUE)

  ageclassdat_L <- ageclass_call(age_formula = ageclass_formula,
                                 obs_aux = format_cjs$l_aux_df,
                                 # max_a = format_cjs$max_a,
                                 ll = TRUE)

  # marginalize these
  ageclassdat_L$mod_mat_a_beta
  # need to make sure the first five columns are j, s, b, obs_time, ageclass
  # if (sum(colnames(s4t_ch$ch$l_matrix)[1:6] == c("j","s","b","g","obs_time","ageclass")) != 6) {
  #   stop("l_matrix not formed correctly with first five columns: j, s, b, obs_time, ageclass")
  # }

  ## Marginalize l_matrix
  unique_identifier_l <- cbind(ageclassdat_L$mod_mat_a_beta,
                               s4t_ch$ch$l_matrix[,1:6])


  unique_l <- !duplicated(unique_identifier_l)

  # n_L <- rep(1,sum(unique_l))

  suppressMessages(marg_unique_id_l <- unique_identifier_l %>% as.data.frame() %>%
                     dplyr::group_by_all() %>%
                     dplyr::summarize(n = dplyr::n()))



  # redefine ageclassdat_L



  # ageclassdat_L$mod_mat_a_beta <- marg_unique_id_l[,1:ncol(ageclassdat_L$mod_mat_a_beta)] # ageclassdat_L$mod_mat_a_beta[unique_l,]
  # ageclassdat_L$obsageclass <-  marg_unique_id_l$ageclass# ageclassdat_L$obsageclass[unique_l]



  l_matrix_marg <- marg_unique_id_l[,c(1:7) + ncol(ageclassdat_L$mod_mat_a_beta)] # s4t_ch$ch$l_matrix[unique_l,1:5]

  l_matrix_marg <- as.matrix(l_matrix_marg)
  # l_matrix_marg <- cbind(l_matrix_marg,n = n_L)


  # ageclassdat_L <- ageclass_call(age_formula = ageclass_formula,
  #                                obs_aux = marg_unique_id_l,
  #                                max_a = format_cjs$max_a,ll = TRUE)

  ageclassdat_L$mod_mat_a_beta <- as.matrix(marg_unique_id_l[,1:ncol(ageclassdat_L$mod_mat_a_beta)])
  ageclassdat_L$obsageclass <- l_matrix_marg[,6]

  # marginalize m_matrix

  unique_identifier_m <- cbind(ageclassdat_M$mod_mat_a_beta,
                               s4t_ch$ch$m_matrix[,1:8])

  suppressMessages(marg_unique_id_m <- unique_identifier_m %>% as.data.frame() %>%
                     dplyr::group_by_all() %>%
                     dplyr::summarize(n = dplyr::n()))



  m_matrix_marg <- marg_unique_id_m[,c(1:9) + ncol(ageclassdat_M$mod_mat_a_beta)] # s4t_ch$ch$m_matrix[unique_m,1:7]
  # m_matrix_marg <- cbind(m_matrix_marg,n = n_M)
  m_matrix_marg <- as.matrix(m_matrix_marg)

  # ageclassdat_M <- ageclass_call(age_formula = ageclass_formula,
  #                                obs_aux = marg_unique_id_m,
  #                                max_a = format_cjs$max_a,ll = TRUE)

  ageclassdat_M$mod_mat_a_beta <- as.matrix(marg_unique_id_m[,1:ncol(ageclassdat_M$mod_mat_a_beta)])
  ageclassdat_M$obsageclass <- m_matrix_marg[,8]

  # not marginalizing this portion
  ageclass_data <- ageclass_call(age_formula=ageclass_formula,
                                 obs_aux = s4t_ch$ch$obs_aux[stats::complete.cases(s4t_ch$ch$obs_aux),]
                                 #na.omit(s4t_ch$ch$obs_aux),
                                 # max_a = s4t_ch$ch_info$max_a
  )

  ageclassdat_L <- ageclassdat_L
  ageclassdat_M <- ageclassdat_M

  # HERE
  if (fixed_age) {
    ageclass_fit <- fit_ageclass(age_formula = ageclass_formula,s4t_ch = s4t_ch)

    fixed_ageclass_l <- ageclass_nll(par = ageclass_fit$res$par,
                                     max_a = max(s4t_ch$ch_info$set_max_a),
                                     mod_mat_a_beta = ageclassdat_L$mod_mat_a_beta,
                                     ll = FALSE)

    fixed_ageclass_m <- ageclass_nll(par = ageclass_fit$res$par,
                                     max_a = max(s4t_ch$ch_info$set_max_a),,
                                     mod_mat_a_beta = ageclassdat_M$mod_mat_a_beta,
                                     ll = FALSE)


  } else {
    fixed_ageclass_l <- fixed_ageclass_m <- NULL
  }


  # set inits and lower and upper bounds

  if (fixed_age) {
    par <- inits <- c(format_cjs$p_par,format_cjs$theta_par)


    lower <- rep(-Inf,length(inits))

    upper <- rep(Inf,length(inits))


  }else {
    par <- inits <- c(ageclass_data$age_par,format_cjs$p_par,format_cjs$theta_par)

    lower <- ifelse(grepl("a_alpha_",names(inits)) &
                      !grepl("^a_alpha_1$",names(inits)),
                    0,-Inf)

    upper <- rep(Inf,length(inits))


  }


  names(lower) <- names(upper) <- names(inits)

  mod_mat_theta; indices_theta

  res <- stats::optim(par = inits,fn = cdl_alg_marg,method = "L-BFGS-B",
               lower = lower,
               upper = upper,
               m_matrix = m_matrix_marg,
               l_matrix = l_matrix_marg,
               obs_aux = s4t_ch$ch$obs_aux,
               max_t_recap = s4t_ch$ch_info$max_t_recap,
               max_s_rel = s4t_ch$ch_info$max_s_rel,
               sites = s4t_ch$user_defined$sites_config,
               site_path=s4t_ch$ch_info$site_path,
               holdover_config = s4t_ch$user_defined$holdover_config,
               set_min_a = s4t_ch$ch_info$set_min_a,
               set_max_a = s4t_ch$ch_info$set_max_a,
               n_batches  = s4t_ch$ch_info$n_batches,
               batches_list = s4t_ch$ch_info$batches_list,
               mod_mat_theta = format_cjs$mod_mat_theta,
               indices_theta = format_cjs$indices_theta,
               mod_mat_p = format_cjs$mod_mat_p,
               indices_p_obs = format_cjs$indices_p_obs,
               ageclass_data = ageclass_data,
               ageclassdat_L = ageclassdat_L,
               ageclassdat_M = ageclassdat_M,
               fixed_age = fixed_age,
               fixed_ageclass_l = fixed_ageclass_l,
               fixed_ageclass_m = fixed_ageclass_m,
               hessian = TRUE,
               control = list(maxit = maxit,
                              lmm = lmm,
                              ndeps = rep(ndeps,length(inits))))

  if (FALSE) {
    par = inits
    m_matrix = m_matrix_marg
    l_matrix = l_matrix_marg
    obs_aux = s4t_ch$ch$obs_aux
    max_t_recap = s4t_ch$ch_info$max_t_recap
    max_s_rel = s4t_ch$ch_info$max_s_rel
    sites = s4t_ch$user_defined$sites_config
    site_path=s4t_ch$ch_info$site_path
    holdover_config = s4t_ch$user_defined$holdover_config
    set_min_a = s4t_ch$ch_info$set_min_a
    set_max_a = s4t_ch$ch_info$set_max_a
    n_batches  = s4t_ch$ch_info$n_batches
    batches_list = s4t_ch$ch_info$batches_list
    mod_mat_theta = mod_mat_theta
    indices_theta = indices_theta
    mod_mat_p = mod_mat_p
    indices_p_obs = indices_p_obs
    ageclass_data = ageclass_data
    ageclassdat_L = ageclassdat_L
    ageclassdat_M = ageclassdat_M
    fixed_age = fixed_age
    fixed_ageclass_l = fixed_ageclass_l
    fixed_ageclass_m = fixed_ageclass_m
  }


  if (res$convergence != 0) {
    message(paste0("Model did not converge. Code: ",res$convergence,"\n",
                   "Convergence message: ",res$message))
  }


  vc <- tryCatch({
    tmpvc <- solve(res$hessian)
    if (sum(is.nan(tmpvc)) > 0) stop("NaN in vcov")

    tmpvc

  },
  error = function(e) {
    message(paste0("Issue with estimating hessian: ","\n",e))
    message("Returning empty vcov")
    matrix(0,nrow = length(res$par),
           ncol = length(res$par),
           dimnames = list(names(res$par),
                           names(res$par)))



  } # end error for trycatch
  ) # end tryCatch

  estimated_parameters <- data.frame(parameter = names(res$par),
                                     estimate = res$par,
                                     se = sqrt(diag(vc)))
  estimated_parameters$lcl95 <- estimated_parameters$estimate - 1.96 * estimated_parameters$se
  estimated_parameters$ucl95 <- estimated_parameters$estimate + 1.96 * estimated_parameters$se
  estimated_parameters$z_stat <- estimated_parameters$estimate / estimated_parameters$se


  s4t_cjs <- list(estimated_parameters = estimated_parameters,
                  res = res,
                  AIC = res$value + 2 * length(res$par),
                  nll = res$value, k = length(res$par),
                  s4t_ch = s4t_ch,
                  fit = list(p_formula = p_formula,
                             theta_formula = theta_formula,
                             ageclass_formula = ageclass_formula,
                             mod_mat_theta = mod_mat_theta,
                             indices_theta = indices_theta,
                             mod_mat_p = mod_mat_p,
                             indices_p_obs = indices_p_obs,
                             ageclass_data = ageclass_data,
                             fixed_age = NULL
                  ))

  if (fixed_age) {
    s4t_cjs$fit$fixed_age <- list(ageclass_fit=ageclass_fit,
                                  fixed_ageclass_l = fixed_ageclass_l,
                                  fixed_ageclass_m = fixed_ageclass_m)

  }

  class(s4t_cjs) <- "s4t_cjs"

  return(s4t_cjs)
}






#' Fit space-for-time mark-recapture model in Bayesian framework
#'
#' @description
#' Uses Stan to fit the age-specific space-for-time model.
#'
#' @param p_formula an object of class "formula" for the formula for detection probabilities.
#' @param theta_formula an object of class "formula" for the formula for transition probabilities.
#' @param ageclass_formula an object of class "formula" for the ageclass sub-model
#' @param cov_p a `data.frame` or `list` of `data.frame`'s containing the covariates for p. See details.
#' @param cov_theta a `data.frame` or `list` of `data.frame`'s containing the covariates for theta. See details.
#' @param groups a `character` vector containing the names of the covariates that comprise the groups.
#' @param s4t_ch a `s4t_ch` object
#' @param chains an `integer` of the number of chains to run.
#' @param warmup an `integer` of the number of warmup iterations.
#' @param iter an `integer` of the number of total (warmup + actual) iterations.
#' @param fixed_age a `logical` object that determines whether the ageclass model will be run
#'     as a separate model (TRUE) or whether it is estimated along with the CJS model (FALSE).
#' @param ... further arguents to pass to `rstan::sampling`
#'
#' @returns a `s4t_cjs` object.
#' @examples
#' # don't run
#'
#' @export
fit_s4t_cjs_rstan <- function(p_formula,
                                   theta_formula,
                                   ageclass_formula,
                                   cov_p = NULL, cov_theta = NULL,groups = NULL,
                                   s4t_ch,
                                   chains = 3,
                                   warmup = 500,
                                   iter = 1000,
                                   fixed_age = TRUE,
                                   ...) {

  dots <- list(...)

  # fixed_age <- match.arg(fixed_age,choices = c(TRUE,FALSE))

  call <- list(p_formula = p_formula,
               theta_formula = theta_formula,
               cov_p = cov_p,
               cov_theta = cov_theta,
               groups =groups,
               s4t_ch = s4t_ch,
               chains = chains,
               warmup = warmup,
               iter = iter,
               fixed_age = fixed_age,
               dots = dots)

  format_cjs <- format_s4t_cjs(p_formula =p_formula, theta_formula = theta_formula,
                               ageclass_formula = ageclass_formula,
                               cov_p = cov_p,cov_theta = cov_theta,s4t_ch = s4t_ch)


  max_t_recap <- format_cjs$max_t_recap
  set_min_a <- format_cjs$set_min_a
  set_max_a <- format_cjs$set_max_a
  max_a_overall <- max(set_max_a)
  N_groups <- format_cjs$N_groups
  tmp_obs_aux <-format_cjs$obs_aux[!is.na(format_cjs$obs_aux[,"ageclass"]),]
  ageclassdat <- ageclass_call(age_formula = ageclass_formula,
                               obs_aux = tmp_obs_aux
                               # max_a = format_cjs$max_a,ll = TRUE
  )
  ageclassdat_M <- ageclass_call(age_formula = ageclass_formula,
                                 obs_aux = format_cjs$m_aux_df,
                                 # max_a = format_cjs$max_a,
                                 ll = TRUE)

  ageclassdat_L <- ageclass_call(age_formula = ageclass_formula,
                                 obs_aux = format_cjs$l_aux_df,
                                 # max_a = format_cjs$max_a,
                                 ll = TRUE)

  min_ageclass_mat <- s4t_ch$ch_info$min_ageclass_mat
  max_ageclass_mat <- s4t_ch$ch_info$max_ageclass_mat
  init_rel_j <- s4t_ch$ch_info$init_rel_j
  init_rel_times <- s4t_ch$ch_info$init_rel_times

  ## Marginalize l_matrix
  unique_identifier_l <- cbind(ageclassdat_L$mod_mat_a_beta,
                               s4t_ch$ch$l_matrix[,1:6])


  unique_l <- !duplicated(unique_identifier_l)

  # n_L <- rep(1,sum(unique_l))

  suppressMessages(marg_unique_id_l <- unique_identifier_l %>% as.data.frame() %>%
                     dplyr::group_by_all() %>%
                     dplyr::summarize(n = dplyr::n()))



  # redefine ageclassdat_L



  # ageclassdat_L$mod_mat_a_beta <- marg_unique_id_l[,1:ncol(ageclassdat_L$mod_mat_a_beta)] # ageclassdat_L$mod_mat_a_beta[unique_l,]
  # ageclassdat_L$obsageclass <-  marg_unique_id_l$ageclass# ageclassdat_L$obsageclass[unique_l]



  l_matrix_marg <- marg_unique_id_l[,c(1:7) + ncol(ageclassdat_L$mod_mat_a_beta)] # s4t_ch$ch$l_matrix[unique_l,1:5]

  l_matrix_marg <- as.matrix(l_matrix_marg)
  # l_matrix_marg <- cbind(l_matrix_marg,n = n_L)


  # ageclassdat_L <- ageclass_call(age_formula = ageclass_formula,
  #                                obs_aux = marg_unique_id_l,
  #                                max_a = format_cjs$max_a,ll = TRUE)

  ageclassdat_L$mod_mat_a_beta <- as.matrix(marg_unique_id_l[,1:ncol(ageclassdat_L$mod_mat_a_beta)])
  ageclassdat_L$obsageclass <- l_matrix_marg[,6]

  # marginalize m_matrix

  unique_identifier_m <- cbind(ageclassdat_M$mod_mat_a_beta,
                               s4t_ch$ch$m_matrix[,1:8])

  suppressMessages(marg_unique_id_m <- unique_identifier_m %>% as.data.frame() %>%
                     dplyr::group_by_all() %>%
                     dplyr::summarize(n = dplyr::n()))



  m_matrix_marg <- marg_unique_id_m[,c(1:9) + ncol(ageclassdat_M$mod_mat_a_beta)] # s4t_ch$ch$m_matrix[unique_m,1:7]

  m_matrix_marg <- as.matrix(m_matrix_marg)


  ageclassdat_M$mod_mat_a_beta <- as.matrix(marg_unique_id_m[,1:ncol(ageclassdat_M$mod_mat_a_beta)])
  ageclassdat_M$obsageclass <- m_matrix_marg[,8]

  if (fixed_age) {
    # currently just the ML fit
    ageclass_fit <- fit_ageclass(age_formula = ageclass_formula,s4t_ch = s4t_ch)

    fixed_ageclass_l <- ageclass_nll(par = ageclass_fit$res$par,
                                     max_a = max(set_max_a),
                                     mod_mat_a_beta = ageclassdat_L$mod_mat_a_beta,
                                     ll = FALSE)

    fixed_ageclass_m <- ageclass_nll(par = ageclass_fit$res$par,
                                     max_a = max(set_max_a),
                                     mod_mat_a_beta = ageclassdat_M$mod_mat_a_beta,
                                     ll = FALSE)


  } else {
    fixed_ageclass_l <- matrix(0,
                               nrow = nrow(ageclassdat_L$mod_mat_a_beta),
                               ncol = max(set_max_a))
    fixed_ageclass_m <- matrix(0,
                               nrow = nrow(ageclassdat_M$mod_mat_a_beta),
                               ncol = max(set_max_a))
  }



  tmp <- max(unlist(lapply(s4t_ch$ch_info$site_path,FUN = length)))
  site_path_mat <- matrix(0,nrow = length(s4t_ch$ch_info$site_path),ncol = tmp)
  for (i in 1:length(s4t_ch$ch_info$site_path)) {
    site_path_mat[i,] = c(s4t_ch$ch_info$site_path[[i]],rep(0,tmp - length(s4t_ch$ch_info$site_path[[i]])))
  }

  tmp <- max(unlist(lapply(s4t_ch$ch_info$batches_list,FUN = length)))
  batches_list_mat <- matrix(0,nrow = length(s4t_ch$ch_info$batches_list),ncol = tmp)
  for (i in 1:length(s4t_ch$ch_info$batches_list)) {
    batches_list_mat[i,] = c(s4t_ch$ch_info$batches_list[[i]],rep(0,tmp - length(s4t_ch$ch_info$batches_list[[i]])))
  }

  #
  # tmp <- max(unlist(lapply(s4t_ch$ch_info$init_rel_times,FUN = length)))
  # init_rel_times_mat <- matrix(0,
  #                              ncol = max(init_rel_j),
  #                              nrow = length(s4t_ch$ch_info$init_rel_times))
  # for (j in init_rel_j) {
  #   init_rel_times_mat[,j] = c(init_rel_times[[j]],rep(0,tmp - length(init_rel_times[[j]])))
  # }


  inits_theta = array(0,dim = c(g = N_groups,
                                j = s4t_ch$ch_info$n_stations,
                                s = max(s4t_ch$ch_info$max_t),
                                b = s4t_ch$ch_info$n_batches,
                                a1 = max_a_overall,
                                a2 = max_a_overall))


  inits_Theta = array(0,dim = c(g = N_groups,
                                j = s4t_ch$ch_info$n_stations,
                                s = max(s4t_ch$ch_info$max_t),
                                b = s4t_ch$ch_info$n_batches,
                                a1 = max_a_overall,
                                a2 = max_a_overall + 1))



  inits_lambda_array = array(0,dim = c(g = N_groups,
                                       j = s4t_ch$ch_info$n_stations,
                                       k = s4t_ch$ch_info$n_stations,
                                       t = max(s4t_ch$ch_info$max_t),
                                       b = s4t_ch$ch_info$n_batches,
                                       s = max(s4t_ch$ch_info$max_t),
                                       a1 = max_a_overall))

  inits_p_obs = array(0,dim = c(g = N_groups,
                                k = s4t_ch$ch_info$n_stations,
                                b = s4t_ch$ch_info$n_batches,
                                t = max(format_cjs$indices_p_obs[,"t"]), # check
                                a1 = max(set_max_a),
                                a2 = max(set_max_a))
  )

  inits_p_obs[,s4t_ch$ch_info$last_sites,,,,] <- 1


  inits_chi_array <- array(0,dim = c(g = N_groups,
                                     j = s4t_ch$ch_info$n_stations,
                                     b = s4t_ch$ch_info$n_batches,
                                     s = max(s4t_ch$ch_info$max_t),
                                     a1 = max_a_overall))

  inits_chi_array[,s4t_ch$ch_info$last_sites,,,] <- 1


  overall_surv <- data.frame(j = as.integer(),k = as.integer(),
                             a1 = as.integer(),s= as.integer(),
                             b = as.integer(),
                             g = as.integer())
  not_last_sites = c(1:s4t_ch$ch_info$n_stations)[-s4t_ch$ch_info$last_sites]
  next_site = apply(s4t_ch$user_defined$sites_config,
                    MARGIN = 1,FUN = function(x) ifelse(is.na(which(x == 1)[1]),0,which(x == 1)[1]))


  # overall survival
  for (g in 1:N_groups) {
    for (j in not_last_sites) {
      k = next_site[j];


      for (s in 1:format_cjs$max_s_rel[j]) {
        tmp_min_a <- min_ageclass_mat[j,s]
        tmp_max_a <- max_ageclass_mat[k,s]

        if (is.na(tmp_min_a) | is.na(tmp_max_a)) {
          next()
          # if it is NA, then skip, because it is likely a
          # site-time combination at a release site when no individuals were released
        }


        for (a1 in tmp_min_a:tmp_max_a) {
          for (b in format_cjs$batches_list[[j]]) {
            overall_surv <- rbind(overall_surv,c(j=j,k=k,
                                                 a1=a1,s=s,
                                                 b=b,g=g))

          }

        }
      }

    }
  } # end

  colnames(overall_surv) <- c("j","k","a1","s","b","g")


  cohort_surv <- data.frame(a1 = as.integer(),
                            a2 = as.integer(),
                            s = as.integer(),t = as.integer(),
                            j = as.integer(),k = as.integer(),
                            b = as.integer(),
                            g = as.integer())

  # cohort survival
  for (g in 1:N_groups) {
    for (j in not_last_sites) {
      k = next_site[j];
      for (b in format_cjs$batches_list[[j]]) {
        for (s in 1:format_cjs$max_s_rel[j]) {
          tmp_min_a <- min_ageclass_mat[j,s]
          tmp_max_a <- max_ageclass_mat[k,s]

          if (is.na(tmp_min_a) | is.na(tmp_max_a)) {
            next()
            # if it is NA, then skip, because it is likely a
            # site-time combination at a release site when no individuals were released
          }


          for (a1 in tmp_min_a:tmp_max_a) {
            tmp_max_t <- min(c(max_t_recap[k],
                               s + (set_max_a[k] - a1)))

            for (t in (s):(tmp_max_t)) {
              a2 <- a1 + t - s

              tmp_upperage <- min(c(t - s + a1,tmp_max_a)); tmp_upperage

              # cohort_surv[paste0("cohort_surv[",a1,",",a2,",",s,",",t,",",j,",",k,",",b,"]")] <-
              #   qlogis(Theta[a1,a2,s,j,k,b])
              #
              cohort_surv <- rbind(cohort_surv,c(a1=a1,
                                                 a2=a2,
                                                 s=s,
                                                 t=t,
                                                 j=j,
                                                 k=k,
                                                 b=b,
                                                 g=g))

            } # t
          } # s
        } # a1
      } #  b
    } #j
  }
  colnames(cohort_surv) <- c("a1","a2","s","t","j","k","b","g")



  N_overall_surv <- nrow(overall_surv)
  N_cohort_surv <- nrow(cohort_surv)
  # print(broodyear_surv); print(N_broodyear_surv)

  max_a_overall <- max(s4t_ch$ch_info$set_max_a)
  input_data <- list(max_a_overall = max_a_overall,
                     set_min_a = s4t_ch$ch_info$set_min_a,
                     set_max_a = s4t_ch$ch_info$set_max_a,
                     mod_mat_a_r = nrow(ageclassdat$mod_mat_a_beta),
                     mod_mat_a_c = ncol(ageclassdat$mod_mat_a_beta),
                     N_a_parbeta = ncol(ageclassdat$mod_mat_a_beta)-1,
                     N_obsageclass = length(ageclassdat$obsageclass),
                     obsageclass = ageclassdat$obsageclass,
                     mod_mat_a_beta = as.matrix(ageclassdat$mod_mat_a_beta[,-1]),

                     N_groups = N_groups,
                     N_m = nrow(m_matrix_marg),
                     N_l = nrow(l_matrix_marg),
                     N_j = ncol(s4t_ch$user_defined$sites_config),
                     N_k = length(s4t_ch$ch_info$recap_sites),
                     max_t = max(s4t_ch$ch_info$max_t),#max(c(s4t_ch$ch_info$max_s_rel,s4t_ch$ch_info$max_t_recap)),
                     N_stations = s4t_ch$ch_info$n_stations,
                     N_batches = s4t_ch$ch_info$n_batches,
                     N_recap_sites = length(s4t_ch$ch_info$recap_sites),
                     N_last_sites = length(s4t_ch$ch_info$last_sites),
                     # N_recap_sites_not_last = length(s4t_ch$ch_info$recap_sites_not_last),
                     N_not_last_sites = length(c(1:s4t_ch$ch_info$n_stations)[-s4t_ch$ch_info$last_sites]),

                     max_t_p = max(format_cjs$indices_p_obs[,"t"]), # CHECK
                     N_site_path_length3 = length(which(unlist(lapply(format_cjs$site_path,length)) >= 3)),
                     N_not_last_sites_rev = length(c(1:s4t_ch$ch_info$n_stations)[-s4t_ch$ch_info$last_sites]),
                     max_s_rel = s4t_ch$ch_info$max_s_rel,
                     max_t_recap = s4t_ch$ch_info$max_t_recap,
                     next_site = apply(s4t_ch$user_defined$sites_config,
                                       MARGIN = 1,FUN = function(x) ifelse(is.na(which(x == 1)[1]),0,which(x == 1)[1])),
                     N_knownage_m = length(which(!is.na(m_matrix_marg[,8]))),
                     N_knownage_l = length(which(!is.na(l_matrix_marg[,6]))),
                     knownage_m = which(!is.na(m_matrix_marg[,8])),
                     knownage_l = which(!is.na(l_matrix_marg[,6])),
                     N_unknownage_l = length(which(is.na(l_matrix_marg[,6]))),
                     N_unknownage_m = length(which(is.na(m_matrix_marg[,8]))),
                     unknownage_l = which(is.na(l_matrix_marg[,6])),
                     unknownage_m = which(is.na(m_matrix_marg[,8])),
                     recap_sites = s4t_ch$ch_info$recap_sites,
                     last_sites = as.matrix(s4t_ch$ch_info$last_sites),
                     # recap_sites_not_last = s4t_ch$ch_info$recap_sites_not_last,
                     not_last_sites = c(1:s4t_ch$ch_info$n_stations)[-s4t_ch$ch_info$last_sites],
                     site_path_length3 = as.matrix(which(unlist(lapply(format_cjs$site_path,length)) >= 3)),
                     not_last_sites_rev = rev(c(1:s4t_ch$ch_info$n_stations)[-s4t_ch$ch_info$last_sites]),
                     batches_list_len = unlist(lapply(s4t_ch$ch_info$batches_list,length)),
                     site_path_len = unlist(lapply(s4t_ch$ch_info$site_path,length)),
                     N_overall_surv = N_overall_surv,
                     N_cohort_surv = N_cohort_surv,
                     m_matrix = m_matrix_marg[,c(1:7,9)],
                     l_matrix = l_matrix_marg[,c(1:5,7)],

                     N_theta_r = nrow(format_cjs$mod_mat_theta),
                     N_theta_c = ncol(format_cjs$mod_mat_theta),
                     N_p_r = nrow(format_cjs$mod_mat_p),
                     N_p_c = ncol(format_cjs$mod_mat_p),
                     N_theta_indices_r = nrow(format_cjs$indices_theta),
                     N_theta_indices_c = ncol(format_cjs$indices_theta),
                     N_cohort_surv_r = nrow(cohort_surv),
                     N_cohort_surv_c = ncol(cohort_surv),
                     N_overall_surv_r = nrow(overall_surv),
                     N_overall_surv_c = ncol(overall_surv),
                     N_theta_par = length(format_cjs$theta_par),
                     N_p_indices_r = nrow(format_cjs$indices_p_obs),
                     N_p_indices_c = ncol(format_cjs$indices_p_obs),
                     N_p_par = length(format_cjs$p_par),

                     mod_mat_theta = format_cjs$mod_mat_theta,
                     mod_mat_p = format_cjs$mod_mat_p,
                     indices_theta = format_cjs$indices_theta,
                     indices_p_obs = format_cjs$indices_p_obs,
                     indices_overall_surv = overall_surv,
                     indices_cohort_surv = cohort_surv,

                     batches_list = batches_list_mat,
                     site_path = site_path_mat,

                     mod_mat_a_L_r = nrow(ageclassdat_L$mod_mat_a_beta),
                     mod_mat_a_L_c = ncol(ageclassdat_L$mod_mat_a_beta),
                     mod_mat_a_M_r = nrow(ageclassdat_M$mod_mat_a_beta),
                     mod_mat_a_M_c = ncol(ageclassdat_M$mod_mat_a_beta),

                     N_obsageclass_L = nrow(ageclassdat_L$mod_mat_a_beta),
                     N_obsageclass_M = nrow(ageclassdat_M$mod_mat_a_beta),
                     mod_mat_a_L = as.matrix(ageclassdat_L$mod_mat_a_beta[,-1]),
                     mod_mat_a_M = as.matrix(ageclassdat_M$mod_mat_a_beta[,-1]),
                     obsageclass_L = ifelse(is.na(l_matrix_marg[,6]),0,l_matrix_marg[,6]),
                     obsageclass_M = ifelse(is.na(m_matrix_marg[,8]),0,m_matrix_marg[,8]),

                     inits_theta = inits_theta,
                     inits_Theta = inits_Theta,
                     inits_p_obs = inits_p_obs,
                     inits_lambda_array = inits_lambda_array,
                     inits_chi_array = inits_chi_array,
                     min_ageclass_mat = ifelse(is.na(min_ageclass_mat),
                                               0,
                                               min_ageclass_mat),
                     max_ageclass_mat = ifelse(is.na(max_ageclass_mat),
                                               0,
                                               max_ageclass_mat)
  )

  if (fixed_age) {
    # RUN RES
    input_data[["fixed_ageclass_l"]] <- fixed_ageclass_l
    input_data[["fixed_ageclass_m"]] <- fixed_ageclass_m
    # res <- rstan::stan(file='stan_version/stan_model_s4t_cjs_marg_fixedage_draft6d_04_16_2025.stan', data=input_data,chains = chains,
    #             warmup = warmup,
    #             iter = iter,
    #             ...)
    res <- rstan::sampling(stanmodels$s4t_cjs_fixedage_draft6d,
                           data=input_data,chains = chains,
                           warmup = warmup,
                           iter = iter,
                           ...)


  } else {

    res <- rstan::sampling(stanmodels$s4t_cjs_draft6d,
                           data=input_data,chains = chains,
                           warmup = warmup,
                           iter = iter,
                           ...)

  }


  sum_res <- rstan::summary(res)

  estimated_parameters <- sum_res$summary


  overall_ests <- estimated_parameters[grepl("overall_surv",rownames(estimated_parameters)),]

  colnames(overall_surv) <- c("j","k","a1","t","b","g")

  overall_surv <- cbind(overall_surv,overall_ests)

  cohort_ests <- estimated_parameters[grepl("cohort_surv",rownames(estimated_parameters)),]

  colnames(cohort_surv) <- c("a1","a2","s","t","j","k","b","g")

  cohort_surv <- cbind(cohort_surv,cohort_ests)

  s4t_cjs_rstan <- list(estimated_parameters = estimated_parameters,
                        overall_surv = overall_surv,
                        cohort_surv = cohort_surv,
                        res = res,
                        call = call,
                        fit = list(p_formula = p_formula,
                                   theta_formula = theta_formula,
                                   ageclass_formula = ageclass_formula,
                                   input_data = input_data,
                                   fixed_age = NULL
                        ))

  if (fixed_age) {
    s4t_cjs_rstan$fit$fixed_age <- list(ageclass_fit,fixed_ageclass_l,fixed_ageclass_m)
  }

  class(s4t_cjs_rstan) <- "s4t_cjs_rstan"

  return(s4t_cjs_rstan)

}


