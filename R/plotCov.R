

#' Effect plot for s4t_cjs_rstan object
#'
#' add description
#'
#'
#' @param mod a `s4t_cjs_rstan` object
#' @param par_formula a `character` of the parameter
#' @param conditional NOTE USED
#' @param ... Filters
#' @returns a ggplot2 type figure of the transition probabilities
plotCov <- function(mod, par_formula, conditional = TRUE, ...) {
  if (!is(mod,"s4t_cjs_rstan")) stop("mod must be class s4t_cjs_ml or s4t_cjs_rstan")

  if (param %in% c("a1","a2","s","t","j","k","r","g")) stop("use plotTheta or plotSurvival for this param")

  # param # check that it is a character of length 1
  # par_formula <- ~ temp

# mod$call$p_formula

  format_theta_partial <- format_theta_formula(theta_formula = par_formula,
                               groups = mod$call$groups,
                               s4t_ch = mod$fit$s4t_ch)

  format_theta_full <- format_theta_formula(theta_formula = mod$call$theta_formula,
                                               groups = mod$call$groups,
                                               s4t_ch = mod$fit$s4t_ch)



  full_mod_mat <- format_theta_full$mod_mat_theta
  full_mod_mat[,2:ncol(full_mod_mat)] <- 0

  update_cols <- colnames(format_theta_partial$mod_mat_theta)
  if (length(setdiff(colnames(format_theta_partial$mod_mat_theta),colnames(format_theta_full$mod_mat_theta))) > 0) {
    non_matching_cols <- setdiff(colnames(format_theta_partial$mod_mat_theta),colnames(format_theta_full$mod_mat_theta))

    for (i in 1:length(non_matching_cols)) {
      non_matching_cols[i]
      tmp_vars <- stringr::str_split(non_matching_cols[i],":")[[1]]
      # (?=.*red)(?=.*car)
      match_vars <- colnames(format_theta_full$mod_mat_theta)[grep(paste0("(?=.*",tmp_vars,")",collapse = ""),colnames(format_theta_full$mod_mat_theta),perl = TRUE)]


      match_vars <- match_vars[grepl("^:*$",gsub(paste0(tmp_vars,collapse="|"),"",match_vars))]

      if (length(match_vars) > 1) stop("matching more than one columns in mod_mat_theta")

      update_cols[update_cols==non_matching_cols[i]] <- match_vars

    }
  }

  full_mod_mat[,update_cols] <- format_theta_partial$mod_mat_theta


  # tmp_theta <- mod$estimated_parameters[grepl("^theta_",rownames(mod$estimated_parameters)),c("mean","sd")]



  # estimate covariance matrix (to approximate CIs)
  # if (FALSE) { # approximate
  #   theta_parnames <- mod$original_units$compare_parnames[grepl("theta_params",mod$original_units$compare_parnames)]
  #   tmp_posterior <- do.call(cbind,rstan::extract(mod$res,theta_parnames))
  #   cov_mat <- cov(tmp_posterior)
  #
  #   cov_theta_est <- full_mod_mat %*% tmp_theta[,"mean"]
  #   cov_theta_se <- diag(full_mod_mat %*% cov_mat %*% t(full_mod_mat))
  #
  # }


  ## add option for this to be less accurate and faster
  tmp_Theta_preds <- predictTheta.s4t_cjs_rstan(mod,newdata = full_mod_mat)

  if (ncol(format_theta_full$tmp_indices_theta)>8) {
    ncols_to_add <- ncol(format_theta_full$tmp_indices_theta) - 8
    tmp_Theta_preds <- cbind(tmp_Theta_preds,format_theta_full$tmp_indices_theta[,9:ncol(format_theta_full$tmp_indices_theta)])
    colnames(tmp_Theta_preds)[(1 + ncol(tmp_Theta_preds)-ncols_to_add):ncol(tmp_Theta_preds)] <- colnames(format_theta_full$tmp_indices_theta)[9:(8 + ncols_to_add)]
  }

  # rename the first element in the formula
  # colnames(tmp_Theta_preds)[colnames(tmp_Theta_preds) == as.character(par_formula[[2]])] <- "element1"

  tmp_elements <- (as.character(par_formula[[2]]))
  tmp_elements <- tmp_elements[grepl("[A-Za-z0-9]",tmp_elements)]

  # formulas are apparently not that easy to deal with
  if (sum(grepl("[*]|:",tmp_elements)) > 0) {
    new_tmp <- unlist(stringr::str_split(tmp_elements,pattern = "[*]|:"))
    new_tmp <- gsub(" ","",new_tmp)
    tmp_elements <- new_tmp
  }

  num_elements <- length(tmp_elements)

  # tmp_Theta_preds <- cbind(tmp_Theta_preds,"element1" = tmp_Theta_preds[,as.character(par_formula[[2]])])

  # num_elements <- length(par_formula) - 1
  for (i in 1:num_elements) {
    tmp_Theta_preds <- cbind(tmp_Theta_preds, "placeholder"= tmp_Theta_preds[,tmp_elements[i]])
    colnames(tmp_Theta_preds)[colnames(tmp_Theta_preds) == "placeholder"] <- paste0("element",i)

  }


  if (num_elements == 1) {
    p <- tmp_Theta_preds %>%
      # dplyr::filter(...) %>%
      dplyr::filter(j == 1) %>%
      ggplot2::ggplot(ggplot2::aes(element1,cohort_surv_mean, color = factor(a1))) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = cohort_surv_lcl,ymax = cohort_surv_ucl),width = 0.1)
  } else if (num_elements == 2) {
    p <- tmp_Theta_preds %>%
      dplyr::filter(...) %>%
      # dplyr::filter(j == 1) %>%
      ggplot2::ggplot(ggplot2::aes(element1,cohort_surv_mean)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = cohort_surv_lcl,ymax = cohort_surv_ucl),width = 0.1) +
      ggplot2::facet_wrap(~element2)
  }else if (num_elements == 3) {
    p <- tmp_Theta_preds %>%
      dplyr::filter(...) %>%
      # dplyr::filter(j == 1) %>%
      ggplot2::ggplot(ggplot2::aes(element1,cohort_surv_mean, color = factor(element2))) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(aes(ymin = cohort_surv_lcl,ymax = cohort_surv_ucl),width = 0.1) +
      ggplot2::facet_wrap(~element3)
  }

  # plot(p)
  return(p)


  # if (mod$call$theta_formula)

  # indices_theta_original <- mod$original_units$indices_theta_original
  # indices_theta_incCov <- indices_theta_original[,1:(ncol(indices_theta_original) - 8)] %>%
  #   dplyr::mutate(a1 = factor(a1),
  #                 a2 = factor(a2),
  #                 s = factor(s),
  #                 t = factor(t),
  #                 j = factor(j),
  #                 k = factor(k),
  #                 r = factor(r),
  #                 g = factor(g))

#
#   mod_mat_theta <- mod$fit$input_data$mod_mat_theta
#
#
#   if (conditional) {
#     # if "by" is empty, then plot over the average (all cov should be centered and scaled):
#     if (is.null(by)) {
#       # cov_avg <- mod$fit$input_data
#       tmp_mod_mat <- mod_mat_theta[1:length(unique(mod_mat_theta[,param])),]
#       tmp_mod_mat[,2:ncol(tmp_mod_mat)] <- 0
#       tmp_mod_mat[,param] <- unique(mod_mat_theta[,param])
#
#
#       tmp_theta <- mod$estimated_parameters[grepl("^theta_",rownames(mod$estimated_parameters)),c("mean","sd")]
#
#
#       theta_parnames <- mod$original_units$compare_parnames[grepl("theta_params",mod$original_units$compare_parnames)]
#       tmp_posterior <- do.call(cbind,rstan::extract(mod$res,theta_parnames))
#       cov_mat <- stats::cov(tmp_posterior)
#
#       cov_theta_est <- tmp_mod_mat %*% tmp_theta[,"mean"]
#       cov_theta_se <- diag(tmp_mod_mat %*% cov_mat %*% t(tmp_mod_mat))
#
#       cov_theta_df <- data.frame(estimate = cov_theta_est,
#                                  se = cov_theta_se,
#                                  lcl95 = cov_theta_est - 1.96 * cov_theta_se,
#                                  ucl95 = cov_theta_est + 1.96 * cov_theta_se)
#       cov_theta_df <- cbind(cov_theta_df, param = unique(mod_mat_theta[,param]))
#
#       cov_theta_df %>%
#         ggplot2::ggplot(aes(param, estimate)) +
#         ggplot2::geom_point() + ggplot2::geom_errorbar(aes(ymin = lcl95, ymax = ucl95),width = 0.2)
#
#     }
#     # if "by" is empty, then plot over the average (all cov should be centered and scaled):
#     if (is.null(by)) {
#       # cov_avg <- mod$fit$input_data
#       tmp_mod_mat <- mod_mat_theta[1:length(unique(mod_mat_theta[,param])),]
#       tmp_mod_mat[,2:ncol(tmp_mod_mat)] <- 0
#       tmp_mod_mat[,param] <- unique(mod_mat_theta[,param])
#
#
#       tmp_theta <- mod$estimated_parameters[grepl("^theta_",rownames(mod$estimated_parameters)),c("mean","sd")]
#
#
#       theta_parnames <- mod$original_units$compare_parnames[grepl("theta_params",mod$original_units$compare_parnames)]
#       tmp_posterior <- do.call(cbind,rstan::extract(mod$res,theta_parnames))
#       cov_mat <- stats::cov(tmp_posterior)
#
#       cov_theta_est <- tmp_mod_mat %*% tmp_theta[,"mean"]
#       cov_theta_se <- diag(tmp_mod_mat %*% cov_mat %*% t(tmp_mod_mat))
#
#       cov_theta_df <- data.frame(estimate = cov_theta_est,
#                                  se = cov_theta_se,
#                                  lcl95 = cov_theta_est - 1.96 * cov_theta_se,
#                                  ucl95 = cov_theta_est + 1.96 * cov_theta_se)
#       cov_theta_df <- cbind(cov_theta_df, param = unique(mod_mat_theta[,param]))
#
#       cov_theta_df %>%
#         ggplot2::ggplot(aes(param, estimate)) +
#         ggplot2::geom_point() + ggplot2::geom_errorbar(aes(ymin = lcl95, ymax = ucl95),width = 0.2)
#
#     } else {
#       tmp_mod_mat <- mod_mat_theta %>%
#         as.data.frame() %>%
#         dplyr::select(`(Intercept)`,
#                       matches(paste0("^",param,"$")),
#                       matches(paste0("^",by,"$",collapse = "|")))
#       tmp_mod_mat[,2:ncol(tmp_mod_mat)] <- 0
#       tmp_mod_mat[,param] <- unique(mod_mat_theta[,param])
#
#
#     }
#
#   }
#
#   # mod_mat = model_mat_info(mod$call$theta_formula,df = indices_theta_incCov)$mod_mat
#
#   # 1. determine if param is a factor (being used as a factor), use mod_mat to decide
#   # 2. if factor: loop through and manipulate mod_mat so that it turns on each level at a time.
#   #     also need to do something about interactions.
#   # 3. if numeric, loop through and manipulate mod_mat so that it changes the value one step at a time.
#   #     will need to deal with interactions.
#   # 4. make a figure. x-axis will be the covariate/param. The y-axis will be cohort_transition/surv.
#   #     Add facet_wrap for factors and other variables?
#
#   return(NULL)

}


format_theta_formula <- function(theta_formula,
                           groups = NULL,
                           s4t_ch) {



  holdover_config <- s4t_ch$s4t_config$holdover_config
  sites_config <- s4t_ch$s4t_config$sites_config

  recap_sites <- s4t_ch$ch_info$recap_sites #which(colSums(sites) > 0)
  last_sites <- s4t_ch$ch_info$last_sites # unique(unlist(lapply(site_path,FUN = max)))

  n_sites <- s4t_ch$ch_info$n_sites

  n_init_relsite <- s4t_ch$ch_info$n_init_relsite # n_batches
  init_relsite_list <- s4t_ch$ch_info$init_relsite_list # batches_list

  max_s_rel <- as.integer(s4t_ch$ch_info$max_s_rel)
  max_t_recap <- as.integer(s4t_ch$ch_info$max_t_recap)
  set_max_a <- s4t_ch$s4t_config$set_max_a
  set_min_a <- s4t_ch$s4t_config$set_min_a
  first_obs <- s4t_ch$ch_info$first_obs
  first_sites <- s4t_ch$ch_info$first_sites
  site_path <- s4t_ch$ch_info$site_path

  min_ageclass_mat <- s4t_ch$ch_info$min_ageclass_mat
  max_ageclass_mat <- s4t_ch$ch_info$max_ageclass_mat

  ###
  # groups
  if (!is.null(groups)) {
    if (length(setdiff(groups,colnames(s4t_ch$ch$all_aux))) > 0) {
      stop(paste0("groups were not found in obs_aux: ",paste0(setdiff(groups,colnames(s4t_ch$ch$all_aux)),collapse = ", ")))
    }

    # the below needs to go after the covariates are added in.
    if (length(intersect(groups,"")) > 0) stop("group names are also in covariates")

    df_groups <- as.data.frame(as.data.frame(s4t_ch$ch$all_aux)[,groups])

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

    group_names <- apply(distinct_groups_df[,groups],MARGIN = 1,FUN = function(x) paste0(x, collapse = "_"))
  } else {
    N_groups <- 1
    group_id <- as.factor(rep(1,nrow(s4t_ch$ch$all_aux)))
    group_names <- 1
  } # end if statement for groups


  ###


  recap_sites_not_last <- s4t_ch$ch_info$recap_sites_not_last


  # 7 columns because first two are a1 and a2. Third is t. Fourth is j, fifth is k, sixth is batch, seventh is group
  indices_theta <- matrix(NA, nrow = 0,ncol = 8)
  colnames(indices_theta) <- c("a1","a2","s","t","j","k","r","g")


  for (g in 1:N_groups) {
    for (j in 1:(n_sites-1)) {
      for (k in which(sites_config[j,]==1)) { # doesn't need to be a loop

        # sites[,k]
        # making an assumption that transitions between sites are the same
        # for all batches.
        for (r in init_relsite_list[[j]]) {
          if (holdover_config[j,k] == 0) {
            # min_a <- max(c(1,))
            # max_a +
            # max_s_rel[j]
            # s in 1:max_s_rel[j]
            for (s in 1:max_s_rel[j]) {
              t <- s
              # tmp_min_a <- set_min_a[j] # max(1,t - max_s_rel[r] + 1)
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
                                        r = r,
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
                                        r = r,
                                        g = g)

              # drop transitions that exceed max_t_recap[k]
              tmp_recap_t <- tmp_indices$a2 - tmp_indices$a1 + tmp_indices$s
              tmp_keeprows <- tmp_recap_t <= max_t_recap[k]
              tmp_indices <- tmp_indices[tmp_keeprows,]


              indices_theta <- rbind(indices_theta, tmp_indices)
            }
          }
        } # r

      } # k
    } # j
  }




  # mf2 <- model.frame(~ x * z,data = data.frame(x = 1:10,y = 5:14,z = as.character(rep(c(1:3),length=10))))
  # mt2 <- attr(mf2,"terms")
  # x2 <- model.matrix(mt2, mf2, contrasts)
  tmp_indices_theta <- as.data.frame(apply(indices_theta,MARGIN = 2,as.factor))



  ## Convert indices to original units and use that to merge in covariates.

  sites_names <- colnames(s4t_ch$s4t_config$sites_config)
  j_site_df <- data.frame(j = as.character(1:s4t_ch$ch_info$n_sites),site_rel = as.character(sites_names))
  k_site_df <- data.frame(k = as.character(1:s4t_ch$ch_info$n_sites),site_rec = as.character(sites_names))


  init_relsite_df <- data.frame(r = as.character(1:s4t_ch$ch_info$n_sites),init_relsite = as.character(sites_names))

  group_df <- data.frame(g = as.character(1:N_groups),group_name = group_names)

  time_diff <- s4t_ch$ch_info$observed_relative_min_max$obs_min_time - 1
  age_diff <- min(s4t_ch$ch_info$observed_relative_min_max$obs_min_a) - 1


  indices_p_obs_original <- tmp_indices_p_obs %>%
    as.data.frame() %>%
    dplyr::mutate(a1 = as.character(a1),
                  a2 = as.character(a2),
                  s = as.character(s),
                  t = as.character(t),
                  j = as.character(j),
                  k = as.character(k),
                  r = as.character(r),
                  g = as.character(g)) %>%
    dplyr::left_join(j_site_df, by = "j") %>%
    dplyr::left_join(k_site_df, by = "k") %>%
    dplyr::mutate(time_rel = as.integer(s) + time_diff,
                  time_rec = as.integer(t) + time_diff,
                  age_rel = as.integer(a1) + age_diff,
                  age_rec = as.integer(a2) + age_diff) %>%
    dplyr::left_join(init_relsite_df, by = "r") %>%
    # dplyr::left_join(group_df, by = "g") %>%
    dplyr::select(a1 = age_rel,
                  a2 = age_rec,
                  s = time_rel,
                  t = time_rec,
                  j = site_rel,
                  k = site_rec,
                  r = init_relsite,
                  g = g
                  # g = group_name
    ) %>%
    dplyr::mutate(a1 = factor(a1),
                  a2 = factor(a2),
                  s = factor(s),
                  t = factor(t),
                  j = factor(j),
                  k = factor(k),
                  r = factor(r),
                  g = factor(g))

  indices_theta_original <- tmp_indices_theta %>%
    as.data.frame() %>%
    dplyr::mutate(a1 = as.character(a1),
                  a2 = as.character(a2),
                  s = as.character(s),
                  t = as.character(t),
                  j = as.character(j),
                  k = as.character(k),
                  r = as.character(r),
                  g = as.character(g)) %>%
    dplyr::left_join(j_site_df, by = "j") %>%
    dplyr::left_join(k_site_df, by = "k") %>%
    dplyr::mutate(time_rel = as.integer(s) + time_diff,
                  time_rec = as.integer(t) + time_diff,
                  age_rel = as.integer(a1) + age_diff,
                  age_rec = as.integer(a2) + age_diff) %>%
    dplyr::left_join(init_relsite_df, by = "r") %>%
    # dplyr::left_join(group_df, by = "g") %>%
    dplyr::select(a1 = age_rel,
                  a2 = age_rec,
                  s = time_rel,
                  t = time_rec,
                  j = site_rel,
                  k = site_rec,
                  r = init_relsite,
                  g = g # group_name
    ) %>%
    dplyr::mutate(a1 = factor(a1),
                  a2 = factor(a2),
                  s = factor(s),
                  t = factor(t),
                  j = factor(j),
                  k = factor(k),
                  r = factor(r),
                  g = factor(g))



  ## now with thetas

  # cov_theta <- data.frame(t = as.character(1:4),cov = rnorm(4))
  cov_theta <- s4t_ch$cov_df$cov_theta

  if (!is.null(cov_theta)) {
    if (is(cov_theta,"list")) {

      for (i in 1:length(cov_theta)) {

        match_col <- intersect(colnames(cov_theta[[i]]),colnames(indices_theta_original))

        # cov_theta[[i]][,match_col] <- as.factor(cov_theta[[i]][,match_col])

        for (j in match_col) {
          cov_theta[[i]][,j] <- as.factor(cov_theta[[i]][,j])
        }

        indices_theta_original <- dplyr::left_join(indices_theta_original,stats::na.fail(cov_theta[[i]]))
      }


    }
    if (is(cov_theta,"data.frame")) {
      match_col <- intersect(colnames(cov_theta),colnames(indices_theta_original))

      # cov_theta[,match_col] <- as.factor(cov_theta[,match_col])

      for (j in match_col) {
        cov_theta[,j] <- as.factor(cov_theta[,j])
      }


      indices_theta_original <- dplyr::left_join(indices_theta_original,stats::na.fail(cov_theta))

    }

  }

  if (ncol(indices_theta_original) > 8) {
    tmp_indices_theta <- cbind(tmp_indices_theta,indices_theta_original[,9:ncol(indices_theta_original)])
    colnames(tmp_indices_theta) <- colnames(indices_theta_original)
  } else {
    tmp_indices_theta <- tmp_indices_theta
  }

  if (!is.null(groups)) {
    tmp_indices_theta <- dplyr::left_join(tmp_indices_theta,distinct_groups_df, by ="g")
  }






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



  # mod_mat_theta = mod_mat_theta
  #
  # tmp_indices_theta = tmp_indices_theta

  return(list(mod_mat_theta = mod_mat_theta,

              tmp_indices_theta = tmp_indices_theta
              # min_ageclass_mat = min_ageclass_mat,
              # max_ageclass_mat = max_ageclass_mat
  ))
}

tmp_indices_p_obs <- param <- element1 <- cohort_surv_mean <- cohort_surv_lcl <-
  cohort_surv_ucl <- element2 <- lcl95 <- ucl95 <- `(Intercept)` <- matches <- NULL
