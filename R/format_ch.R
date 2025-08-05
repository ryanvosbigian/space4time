
#' @importFrom magrittr %>%



process_ch <- function(obs_ch,obs_aux,removed_sites,sites = NULL) {
  obs_ch <- as.matrix(obs_ch)
  # obs_aux <- as.matrix(obs_aux)
  # obs_ch <- sim.dat$obs_ch
  # obs_lengthyear <- sim.dat$obs_lengthyear
  max_t <- max(obs_ch)

  # (max_t + 1)*(max_t + 1 - 1) * (max_t + 1 - 2)

  tot_entries_m <- sum(apply(obs_ch, MARGIN = 1,FUN = function(x) sum(x != 0) - 1))

  additional_aux_var <- setdiff(colnames(obs_aux),c("id","obs_time","ageclass"))

  m_matrix <- matrix(0,nrow = tot_entries_m, ncol = 8)
  colnames(m_matrix) <- c("j","k","s","t","r","g","obs_time","ageclass")
  l_matrix <- matrix(0,nrow = nrow(obs_ch),ncol = 6)
  colnames(l_matrix) <- c("j","s","r","g","obs_time","ageclass")


  m_aux_df <- as.data.frame(matrix(0,nrow = tot_entries_m, ncol = 2 + length(additional_aux_var)))
  colnames(m_aux_df) <- c("obs_time","ageclass",additional_aux_var)
  l_aux_df <- as.data.frame(matrix(0,nrow = nrow(obs_ch),ncol = 2 + length(additional_aux_var)))
  colnames(l_aux_df) <- c("obs_time","ageclass",additional_aux_var)

  # placeholder, so 1 for now
  m_matrix[,"g"] <- 1
  l_matrix[,"g"] <- 1

  # m_array <- matrix(0,dim = c(i = entries_m, j= ncol(obs_ch),k = ncol(obs_ch), )

  counter <- 1
  for (i in 1:nrow(obs_ch)) {


    obs <- which(obs_ch[i,] != 0)
    not_obs <- which(obs_ch[i,] == 0)

    first_rel <- obs[1]
    last_rel <- obs[length(obs)]

    #
    entries_m <- length(obs) - 1
    entries_l <- 1 # always 1

    l_matrix[i,1:2] <- c(last_rel,obs_ch[i,last_rel])

    # obtain location of first observation, which corresponds to the init_relsite (batch)
    l_matrix[i,3] <- which.min.not.zero(obs_ch[i,])# which(obs_ch[i,]!= 0)[1]

    l_aux_df[i,"obs_time"] <- l_matrix[i,5] <- obs_aux[i,"obs_time"]
    l_aux_df[i,"ageclass"] <- l_matrix[i,6] <- obs_aux[i,"ageclass"]

    if (length(additional_aux_var) > 0) {
      for (k in 1:length(additional_aux_var)) {
        l_aux_df[i,2+k] <- obs_aux[i,additional_aux_var[k]]
      }
    }



    if (entries_m > 0) {
      for (z in 1:entries_m) {
        m_matrix[counter,1:4] <- c(obs[z],obs[z + 1],obs_ch[i,obs[z]],obs_ch[i,obs[z + 1]])
        m_matrix[counter,5] <- which.min.not.zero(obs_ch[i,]) # which(obs_ch[i,]!= 0)[1]

        m_aux_df[counter,"obs_time"] <- m_matrix[counter,7] <- obs_aux[i,"obs_time"]
        m_aux_df[counter,"ageclass"] <- m_matrix[counter,8] <- obs_aux[i,"ageclass"]

        if (length(additional_aux_var) > 0) {
          for (k in 1:length(additional_aux_var)) {
            m_aux_df[counter,2+k] <- obs_aux[i,additional_aux_var[k]]
          }
        }


        counter <- counter + 1
      }
    }
  }

  # only removing removed individuals from the l_matrix, not the m_matrix.
  # need to add warning if individuals are released back into the system.

  # dropping removed individuals (only keeping individuals )
  l_matrix <- l_matrix[removed_sites == 0,]
  l_aux_df <- l_aux_df[removed_sites == 0,]

  return(list(obs_ch,m_matrix, l_matrix,as.data.frame(m_aux_df),as.data.frame(l_aux_df), obs_aux))

}



faster_process_ch <- function(obs_ch,obs_aux,removed_sites,sites = NULL) {
  obs_ch <- as.matrix(obs_ch)
  # obs_aux <- as.matrix(obs_aux)
  # obs_ch <- sim.dat$obs_ch
  # obs_lengthyear <- sim.dat$obs_lengthyear
  max_t <- max(obs_ch)

  # (max_t + 1)*(max_t + 1 - 1) * (max_t + 1 - 2)

  tot_entries_m <- sum(apply(obs_ch, MARGIN = 1,FUN = function(x) sum(x != 0) - 1))

  additional_aux_var <- setdiff(colnames(obs_aux),c("id","obs_time","ageclass"))

  # # m_matrix <- matrix(0,nrow = tot_entries_m, ncol = 8)
  # colnames(m_matrix) <- c("j","k","s","t","r","g","obs_time","ageclass")
  # # l_matrix <- matrix(0,nrow = nrow(obs_ch),ncol = 6)
  # colnames(l_matrix) <- c("j","s","r","g","obs_time","ageclass")


  m_aux_df <- as.data.frame(matrix(0,nrow = tot_entries_m, ncol = 2 + length(additional_aux_var)))
  colnames(m_aux_df) <- c("obs_time","ageclass",additional_aux_var)
  l_aux_df <- as.data.frame(matrix(0,nrow = nrow(obs_ch),ncol = 2 + length(additional_aux_var)))
  colnames(l_aux_df) <- c("obs_time","ageclass",additional_aux_var)

  # # placeholder, so 1 for now
  # m_matrix[,"g"] <- 1
  # l_matrix[,"g"] <- 1

  # m_array <- matrix(0,dim = c(i = entries_m, j= ncol(obs_ch),k = ncol(obs_ch), )


  l_matrix <- apply(obs_ch,MARGIN = 1,FUN = function(ch) {
    last_obs <- sum(ch != 0)
    last_site <- which(ch != 0)[last_obs]
    c(j = last_site,
      s = ch[last_site],
      r = which.min.not.zero(ch),
      g = 1,
      obs_time = NA,
      ageclass = NA)
  }
  )
  l_matrix <- t(l_matrix)
  colnames(l_matrix) <- c("j","s","r","g","obs_time","ageclass")

  l_matrix[,"obs_time"] <- obs_aux$obs_time
  l_matrix[,"ageclass"] <- obs_aux$ageclass

  l_aux_df[,"obs_time"] <- obs_aux$obs_time
  l_aux_df[,"ageclass"] <- obs_aux$ageclass
  l_aux_df[,additional_aux_var] <- obs_aux[,additional_aux_var]



  iterate_m_entries <- apply(obs_ch, MARGIN = 1,FUN = function(x) sum(x != 0) - 1)
  iterate_m_entries_drop0 <- iterate_m_entries[iterate_m_entries > 0]

  recap.num <- unlist(sapply(iterate_m_entries_drop0, FUN = function(x) 1:x))

  iterate_m_ids <- rep(1:nrow(obs_ch),times = iterate_m_entries)
  iterate_m <- cbind(id = iterate_m_ids,
                     recap.num = recap.num)

  tmp_m_matrix <- apply(iterate_m,MARGIN = 1,
        FUN = function(x) {
          ch <- obs_ch[x[1],]

          r <- which.min.not.zero(ch)

          first_site <- which(ch != 0)[1]
          all_site <- which(ch != 0)


          j = all_site[x[2]]
          k = all_site[x[2] + 1]

          matrix(c(j = j,
            k = k,
            s = ch[j],
            t = ch[k],
            r = r,
            g = 1,
            obs_time = NA,
            ageclass = NA,
            id = x[1]),
            nrow = 1,
            ncol = 9)
        }
        )
  tmp_m_matrix <- t(tmp_m_matrix)
  colnames(tmp_m_matrix) <- c("j","k","s","t","r","g","obs_time","ageclass","id")

  tmp_m_matrix[,"obs_time"] <- obs_aux$obs_time[tmp_m_matrix[,"id"]]
  tmp_m_matrix[,"ageclass"] <- obs_aux$ageclass[tmp_m_matrix[,"id"]]

  m_aux_df[,"obs_time"] <- obs_aux$obs_time[tmp_m_matrix[,"id"]]
  m_aux_df[,"ageclass"] <- obs_aux$ageclass[tmp_m_matrix[,"id"]]
  for (j in 1:length(additional_aux_var)) {
    m_aux_df[,additional_aux_var[j]] <- obs_aux[tmp_m_matrix[,"id"],
                                                additional_aux_var[j]]

  }

  m_matrix <- tmp_m_matrix[,c(1:8)]


  colnames(m_matrix) <- c("j","k","s","t","r","g","obs_time","ageclass")





  # only removing removed individuals from the l_matrix, not the m_matrix.
  # need to add warning if individuals are released back into the system.

  # dropping removed individuals (only keeping individuals )
  l_matrix <- l_matrix[removed_sites == 0,]
  l_aux_df <- l_aux_df[removed_sites == 0,]

  return(list(obs_ch,m_matrix, l_matrix,as.data.frame(m_aux_df),as.data.frame(l_aux_df), obs_aux))

}


new_s4t_ch <- function(obs_ch,
                       obs_aux,
                       all_aux,
                       ch_df,
                       set_max_a,
                       s4t_config,
                       observed_relative_min_max,
                       potential_error_log,
                       cov_p = NULL,
                       cov_theta = NULL,
                       call) {

  sites_config <- s4t_config$sites_config
  holdover_config <- s4t_config$holdover_config

  mat_obs_ch <- as.matrix(obs_ch[,colnames(sites_config)])
  removed_sites <- (obs_ch[,ncol(obs_ch)])
  # obs_aux <- as.matrix(obs_aux)


  # message("Starting process_ch")

  ch <- faster_process_ch(obs_ch = mat_obs_ch,obs_aux = obs_aux,removed_sites)

  # message("Finished process_ch")

  m_matrix <- ch[[2]]
  l_matrix <- ch[[3]]

  m_aux_df <- ch[[4]]
  l_aux_df <- ch[[5]]

  site_path <- list()

  for (j in 1:nrow(sites_config)) {

    tmp <- j
    keep_going <- TRUE
    while(keep_going) {
      lastelement <- tmp[length(tmp)]
      if (rowSums(sites_config)[lastelement] == 0) {
        keep_going <- FALSE
      } else {
        tmp <- c(tmp,which(sites_config[lastelement,] == 1))
      }

    }
    site_path[[j]] <- tmp
  }
  site_path


  # # First site in path
  site_1stin_path <- list()
  first_sites <- which(colSums(sites_config) == 0)

  for (j in 1:nrow(sites_config)) {

    f <- j
    keep_going <- TRUE
    while (keep_going) {
      f_int_1st <- intersect(f,first_sites)
      # r_diff_init_relsite <- setdiff(f,init_relsite )
      # intersect(r_diff_init_relsite,first_sites)

      if (length(f_int_1st) == length(f)){

        keep_going <- FALSE
      } else {
        tmp_f1 <- intersect(f,first_sites)
        tmp_f2 <- f
        tmp_f3 <- setdiff(f,first_sites)

        app <- c()
        if (length(tmp_f3) > 0) {
          for (i in 1:length(tmp_f3)) {
            app <-  c(app,which(sites_config[,tmp_f3[i]] == 1))
          }
        }



        f <- c(tmp_f1,app)
      }


    }

    site_1stin_path[[j]] <- f
  }
  site_1stin_path

  # message("checkpoint 1")

  # first_obs <- apply(obs_ch,MARGIN = 1,FUN = function(x) which(x != 0)[1])
  first_obs <- apply(obs_ch[,1:(ncol(obs_ch)-1)],MARGIN = 1,FUN = function(x) which.min.not.zero(x)[1])

  # sort(unique(first_obs))

  # which(sim.dat$obs_ch[,1])


  init_relsite <- which(colSums(sites_config) == 0) # sort(unique(first_obs)) #
  n_init_relsite <- length(init_relsite)# sum(colSums(sites) == 0)

  init_relsite_list <- list()

  first_sites <- which(colSums(sites_config) == 0)


  # message("checkpoint 1.b")

  for (j in 1:nrow(sites_config)) {

    r <- j
    # r <- c()

    keep_going <- TRUE
    while (keep_going) {
      r_int_init_relsite <- intersect(r,init_relsite)
      r_diff_init_relsite <- setdiff(r,init_relsite)

      # NEED to check whether all sites are properly explored

      explored <- TRUE
      for (i in 1:length(r)) {
        # if there are any sites in the path that haven't been explored.
        if (length(setdiff(site_1stin_path[[r[i]]],r)) > 0) {
          explored <- FALSE
        }
      }


      if (length(unique(intersect(r,init_relsite))) == length(unique(r)) &
          explored == TRUE){

        keep_going <- FALSE
      } else {
        # save all init_relsite encountered
        tmp_r1 <- intersect(r,init_relsite)
        tmp_r2 <- r

        # save all the r's that are a first site.
        tmp_r3 <- setdiff(r,first_sites)
        tmp_r4 <- intersect(r,first_sites)

        # loop through all the r's that are not a first site to explore it further
        app_tmp <- c()
        if (length(tmp_r3) > 0) {
          for (i in 1:length(tmp_r3)) {
            app_tmp <-  c(app,which(sites_config[,tmp_r3[i]] == 1))
          }
        }

        r <- unique(sort(c(tmp_r1,tmp_r4,app_tmp))) # unique


      }


    }


    init_relsite_list[[j]] <- r

    # message(paste0("checkpoint 1.c: j = ",j))

  } # end for loop
  init_relsite_list

  # message("checkpoint 2")

  n_sites <- ncol(sites_config)

  # should have as many entries as n_sites
  max_s_rel <- rep(0,n_sites); names(max_s_rel) <- 1:n_sites
  # max_s_rel <- unlist(lapply(split(l_matrix[,"s"],l_matrix[,"j"]),FUN = max)); max_s_rel
  tmp1 <- unlist(lapply(split(l_matrix[,"s"],l_matrix[,"j"]),FUN = max)); tmp1
  tmp2 <- unlist(lapply(split(m_matrix[,"s"],m_matrix[,"j"]),FUN = max)); tmp2
  max_s_rel[names(tmp1)] <- tmp1; max_s_rel

  max_s_rel[names(tmp2)] <- ifelse(max_s_rel[names(tmp2)] > tmp2, max_s_rel[names(tmp2)], tmp2)

  common_occ <- intersect(names(tmp1),names(tmp2))
  if (sum(tmp1[common_occ] != tmp2[common_occ]) > 0) {
    message("Maximum release occasions in l_matrix and m_matrix do not match.")
  }



  # Code to drop observations from the last sites

  # l_matrix <- l_matrix[which(!(l_matrix[,"j"] %in% last_sites)),]

  # max_t for sites with recaps. Typically second to last n_sites (when there is only 1 release site)
  max_t_recap <- rep(0,n_sites); names(max_t_recap) <- 1:n_sites
  tmp <- unlist(lapply(split(m_matrix[,"t"],m_matrix[,"k"]),FUN = max))
  max_t_recap[names(tmp)] <- tmp; max_t_recap


  first_sites
  tmp_nohold <- setdiff(which(colSums(holdover_config)==0),first_sites)

  ### PART TO INVESTIGATE

  ## This takes the max_s_rel and max_t_recap
  ## and makes it so that the max times between sites that individuals
  ## migrate directly (don't holdover) equal to the max time that an individual
  ## was recapped or released

  if (length(tmp_nohold) > 0) {
    # I go over it twice, once in either direction so that
    # the max times will proliferate
    for (h in c(1:length(tmp_nohold),length(tmp_nohold):1)) {
      tmp_nohold[h]
      tmp_prevsites <- which(sites_config[,tmp_nohold[h]] == 1)

      # sapply(c(3,1,3),FUN = function(x) max(c(x,2)))
      max_s_rel[tmp_prevsites] <- sapply(max_s_rel[tmp_prevsites],
                                         FUN = function(x) max(c(x,max_s_rel[tmp_nohold[h]])))

      max_t_recap[tmp_prevsites] <- sapply(max_t_recap[tmp_prevsites],
                                           FUN = function(x) max(c(x,max_t_recap[tmp_nohold[h]])))



    }
  }
  max_s_rel
  max_t_recap


  # message("checkpoint 3")
  ## Determines the minimum and maximum age for each site-time.
  # This is needed, because otherwise, the likelihood will include
  # transitions between ages that are not possible.

  # determine first captures
  first_cap_j  <- apply(obs_ch[,1:(ncol(obs_ch)-1)],MARGIN = 1,
                        FUN = function(x) which(x != 0)[1])
  tmp_first_cap_diff <- setdiff(unique(first_cap_j),init_relsite)
  if (length(tmp_first_cap_diff) > 0) {
    warning(paste0("First observation of some (N=",sum(first_cap_j %in% tmp_first_cap_diff),") individuals are in recap sites: ",
                   paste0(s4t_config$sites_names[tmp_first_cap_diff],collapse = ", ")))
  }
  # warning()

  init_rel_j <- sort(unique(first_cap_j))
  min_init_rel <- rep(0,length = length(init_rel_j)) # is this right??
  max_init_rel <- rep(0,length = length(init_rel_j))

  init_rel_times <- list()

  for (j in init_rel_j) {
    tmp <- obs_ch[first_cap_j == j,j]
    min_init_rel[j] <- min(tmp,na.rm = TRUE)
    max_init_rel[j] <- max(tmp,na.rm = TRUE)

    init_rel_times[[j]] <- sort(unique(tmp)) # obtains release times

  }
  init_rel_j
  min_init_rel
  max_init_rel

  min_ageclass_mat <- max_ageclass_mat <- matrix(data = NA,nrow = nrow(sites_config),
                                                 ncol = max(max_t_recap))

  set_min_a <- observed_relative_min_max$set_min_a
  set_max_a <- observed_relative_min_max$set_max_a

  sites_config
  holdover_config


  # fill out the min-max ageclass matrix for site-time.
  # start with what we know from user-inputs
  for (j in init_rel_j) {
    min_ageclass_mat[j,init_rel_times[[j]]] <- set_min_a[j]
    max_ageclass_mat[j,init_rel_times[[j]]] <- set_max_a[j]
  }

  # next, fill out the min age for all the sites within the times when individuals
  # are released
  for (j in init_rel_j) {
    tmp_next_sites <- setdiff(site_path[[j]],j)
    for (k in tmp_next_sites) {
      min_ageclass_mat[k,init_rel_times[[j]]] <- set_min_a[k]
      max_ageclass_mat[k,init_rel_times[[j]]] <- set_max_a[k]


      for (s in 1:ncol(min_ageclass_mat)) {
        # leave as-is if s is in the initial release times
        if (s %in% init_rel_times[[j]]) next()

        # otherwise, need to determine what is the min and max age an individual
        # can be at the site k during time s
        max_s_rel[k]
        init_rel_times[[j]]

        # obtain the difference in ages
        age_diffs <- s - init_rel_times[[j]]

        # add the difference in ages to the minimum age when released. This gives us
        # the actual difference in ages. Drop (set to NA) values that are less than
        # or equal to zero
        min_age_sk <- ifelse(age_diffs + set_min_a[j] < 1,NA,
                             age_diffs + set_min_a[j])




        if (sum(!is.na(age_diffs)) == 0) stop("Error in determining minimum age for site-times")

        existing_val_min_tmp <- min_ageclass_mat[k,s]

        suppressWarnings(min_ageclass_mat[k,s] <- min(c(min_age_sk,
                                       existing_val_min_tmp),na.rm = TRUE))
        if (is.infinite(min_ageclass_mat[k,s])) min_ageclass_mat[k,s] <- NA

        # max values are simple, just given
        max_ageclass_mat[k,s] <- set_max_a[k]


      }

    }
  }

  # message("checkpoint 4")


  recap_sites <- which(colSums(sites_config) > 0)
  last_sites <- unique(unlist(lapply(site_path,FUN = max)))

  recap_sites_not_last <- recap_sites[!(recap_sites %in% last_sites)]


  not_last_sites <- c(1:n_sites)[-last_sites]


  s4t_ch <- list(ch = list(m_matrix = m_matrix,
                           l_matrix = l_matrix,
                           m_aux_df = m_aux_df,
                           l_aux_df = l_aux_df,
                           all_aux = all_aux),
                 obs_data = list(ch_df = ch_df,
                                 obs_aux = obs_aux,
                                 all_aux = all_aux),
                 call = call,
                 s4t_config = s4t_config,
                 # user_defined = list(sites_config = sites_config,
                 #                     holdover_config = holdover_config,
                 #                     max_a = observed_relative_min_max$orig_max_a,
                 #                     set_max_a = set_max_a,
                 #                     sites_names = rownames(sites_config)),
                 ch_info = list(n_init_relsite = n_init_relsite, # n_batches
                                n_sites = n_sites, # n_stations
                                max_s_rel = max_s_rel,
                                max_t_recap = max_t_recap,
                                min_ageclass_mat = min_ageclass_mat,
                                max_ageclass_mat = max_ageclass_mat,
                                # set_max_a = set_max_a,
                                # set_min_a = observed_relative_min_max$set_min_a,
                                observed_relative_min_max = observed_relative_min_max,

                                first_obs = first_obs,
                                first_sites = first_sites,
                                site_path = site_path,
                                recap_sites = recap_sites,
                                last_sites = last_sites,
                                recap_sites_not_last = recap_sites_not_last,

                                init_relsite_list = init_relsite_list,# batches_list
                                init_rel_j = init_rel_j,
                                init_rel_times = init_rel_times),
                 cov_df = list(cov_p = cov_p,
                               cov_theta = cov_theta),
                 potential_error_log = potential_error_log
  )
  class(s4t_ch) = "s4t_ch"

  # message("checkpoint 5")

  return(s4t_ch)
}


#' Create space-for-time mark-recapture capture history object
#'
#' @description
#' Create space-for-time mark-recapture capture history object.
#'
#' @param ch_df `data.frame` object containing a row for each capture event. See details.
#' @param aux_age_df `data.frame` containing auxiliary data for each individual. See details.
#' @param s4t_config a `s4t_config` object created using `s4t_config()`,
#'     `linear_s4t_config`, or `simplebranch_s4t_config`.
#' @param cov_p a `data.frame` or `list` of `data.frame`'s containing the covariates
#'     for p `a1,a2,j,k,s,t,r,g` indices. See details.
#' @param cov_theta a `data.frame` or `list` of `data.frame`'s containing the covariates
#'     for theta  `a1,a2,j,k,s,t,r,g` indices. See details.
#'
#' @details
#' The capture history data (`ch_df`) must be a `data.frame` (or coercible
#'     to a `data.frame`) with exactly four columns named `id`, `site`,
#'     `time`, and `removed`. Each row is a release and recapture (or observation)
#'     event. The `id` column is the unique identifier for each individual.
#'     The `site` column is the site name, which must correspond to the names
#'     in the `s4t_config` object. The `time` column must either be an integer
#'     for the time period or a Date. If it is a date, it will be converted to
#'     years. The last column is a `logical`(i.e. `TRUE` or `FALSE`) that indicates
#'     whether individuals were removed (i.e. retained) at the event.
#'
#' The auxiliary and age data (`aux_age_df`) must be a `data.frame` (or coercible
#'     to a `data.frame`) that must contain at least three columns named `id`,
#'     `obs_time`, and `ageclass`. Additional columns can be included
#'     that contain data on individuals. The `id` column is the unique identifier
#'     for individuals, `obs_time` is the integer time period (or Date) when the
#'     individual was first observed or when the `ageclass` of the individual
#'     was observed. `ageclass` is the integer age of the individual.
#'     If the `ageclass` was not observed, then `ageclass = NA`, but
#'     `obs_time` must be filled in. `obs_time` should correspond to the time period
#'     of the auxiliary data.
#'
#'
#'
#' Note that individual covariates can be included in the `s4t_cjs_ml` and `s4t_cjs_rstan`
#'     models. These covariates are included in the `aux_age_df` data.
#'
#'
#' @examples
#'
#' ch_df <- data.frame(id = c(1,1,1,
#'                            2,2,
#'                            3,3,
#'                            4,
#'                            5,
#'                            6),
#'                      site = c("A","B","C",
#'                               "A","B",
#'                               "A","C",
#'                               "A",
#'                               "A",
#'                               "A"),
#'                      time = c(1,3,3,
#'                               2,3,
#'                               1,3,
#'                               2,
#'                               1,
#'                               1),
#'                      removed = c(FALSE,FALSE,FALSE,
#'                                  FALSE,FALSE,
#'                                  FALSE,FALSE,
#'                                  FALSE,
#'                                  FALSE,
#'                                  FALSE)
#'                       )
#'
#' aux_age_df <- data.frame(id = 1:6,
#'                           obs_site = rep("A",6),
#'                           ageclass = c(1,2,1,1,2,1),
#'                           obs_time = c(1,2,1,2,1,1),
#'                           Covariate1 = c(3,1,2,1,2,1))
#'
#' site_arr <- linear_s4t_config(sites_names = c("A","B","C"),
#'                               holdover_sites = c("A"),
#'                               min_a = c(1,1,1),
#'                               max_a = c(3,3,3))
#'
#' ch <- s4t_ch(ch_df = ch_df,
#'              aux_age_df = aux_age_df,
#'              s4t_config = site_arr)
#'
#'
#'
#'
#' @export
#'
s4t_ch <- function(ch_df,
                   aux_age_df,
                   s4t_config,
                   cov_p = NULL,
                   cov_theta = NULL) {

  ## Consider changing some of these to messages that result in an error after the whole thing is run

  ch_df <- as.data.frame(ch_df)
  aux_age_df <- as.data.frame(aux_age_df)

  # if (!is(ch_df,"data.frame")) {
  #   stop("ch_df must be a data.frame")
  # }
  #
  # if (!is(aux_age_df,"data.frame")) {
  #   stop("aux_age_df must be a data.frame")
  # }


  if (sum(colnames(ch_df) %in% c("id","site","time","removed")) != 4) {

    stop("ch_df must contain columns id, site, time, and removed")
  }


  if (sum(colnames(aux_age_df) %in% c("obs_time","ageclass","id")) != 3) {

    stop("aux_age_df must contain columns obs_time, ageclass, and id")
  }






  # ch_df$id
  # ch_df$site

  # change type of ch_df$id and ch_df$site to characters
  ch_df$id <- as.character(ch_df$id)
  ch_df$site <- as.character(ch_df$site)




  all_sites_names <- c(s4t_config$sites_names,
                       unlist(s4t_config$sites_to_pool))

  tmp_sites_to_remove <- setdiff(unique(ch_df$site),all_sites_names)

  if (length(tmp_sites_to_remove) > 0) {
    message(paste0("Removing sites from capture history: ",paste(tmp_sites_to_remove,sep = "",collapse = ", ")))

    ch_df <- subset(ch_df,!(ch_df$site %in% tmp_sites_to_remove))
  }

  # ch_df$time

  if (is(ch_df$time,"character")) {
    stop("Convert ch_df$time to an integer or date or datetime class")
  }

  if (is(ch_df$time,"Date") |
      is(ch_df$time,"POSIXct") |
      is(ch_df$time,"POSIXt")) {

    message("Coercing ch_df$time to year")

    ch_df$time <- lubridate::year(ch_df$time)
  }

  # repeat with aux_age_df$obs_time


  # change type of ch_df$id and ch_df$site to characters
  aux_age_df$id <- as.character(aux_age_df$id)
  # aux_age_df$obs_site <- as.character(aux_age_df$obs_site)

  if (is(aux_age_df$obs_time,"character")) {
    stop("Convert aux_age_df$obs_time to an integer or date or datetime class")
  }


  if (is(aux_age_df$obs_time,"Date") |
      is(aux_age_df$obs_time,"POSIXct") |
      is(aux_age_df$obs_time,"POSIXt")) {

    message("Coercing aux_age_df$obs_time to year")

    aux_age_df$obs_time <- lubridate::year(aux_age_df$obs_time)
  }

  ###
  ###
  ###
  ###


  # ch_df$removed

  if (sum(is.na(as.logical(ch_df$removed))) > 0) {
    stop("ch_df$removed not coercible to logical or contains NAs")
  }

  # check that sites identified in ch_df and aux_age_df match what is in sites_names
  unique(ch_df$site)



  aux_age_df <- as.data.frame(aux_age_df)
  # use these to recover time classes and age classes
  obs_min_time <- min_not_zero(ch_df$time) # only use CH data. not aux_age_df$obs_time)
  obs_max_time <- max(c(ch_df$time),na.rm = TRUE) #  only use CH data , not aux_age_df$obs_time



  # the minimum obs age-class will be age-class 1 (even if it is 0, 2, etc.)
  obsaux_min_a <- min(aux_age_df$ageclass,na.rm = TRUE)



  # min_obs_age <- min_not_zero(aux_age_df$ageclass)
  obsaux_max_a <- max(aux_age_df$ageclass,na.rm = TRUE)

  # sets the minimum age_class to 1.
  set_max_a <- s4t_config$set_max_a
  set_min_a <- s4t_config$set_min_a
  obs_max_a <- s4t_config$obs_max_a
  obs_min_a <- s4t_config$obs_min_a

  max_rel_age_class = max(s4t_config$obs_max_a) - obsaux_min_a + 1



  ##
  # age - auxiliary variables

  # aux_age_df
  aux_age_df_abbr <- aux_age_df
  aux_age_df_abbr$ageclass <- aux_age_df_abbr$ageclass - min(obs_min_a) + 1
  aux_age_df_abbr$obs_time <- aux_age_df_abbr$obs_time - obs_min_time + 1


  ## next create obs_ch.

  new_ch_df <- data.frame(id = ch_df$id,
                          site = ch_df$site,
                          time = ch_df$time - obs_min_time + 1,
                          removed = as.logical(ch_df$removed))
  new_ch_df$removed <- as.character(ifelse(new_ch_df$removed,new_ch_df$site,FALSE))


  # pool sites
  sites_to_pool <- s4t_config$sites_to_pool

  # using_sites <- sites_names

  if (!is.null(sites_to_pool)) {
    for (i in 1:length(sites_to_pool)) {

      tmp_sitename <- names(sites_to_pool)[i]



      removingsites_df <- new_ch_df %>%
        as.data.frame() %>%
        dplyr::filter(site %in% c(tmp_sitename,sites_to_pool[[i]])) %>% # sites_to_pool[[i]]
        dplyr::mutate(site = paste0(tmp_sitename)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(id = dplyr::first(id),
                  site = paste0(tmp_sitename),
                  time = min(time),
                  removed = (ifelse(any(removed != "FALSE"),tmp_sitename,"FALSE")))
      # distinct(id,site,time,.keep_all = TRUE)
      # mutate(removed = ifelse(removed != FALSE, paste0(tmp_sitename,"_plus"),removed))

      new_ch_df_dropsites <- new_ch_df %>%
        dplyr::filter(!(site %in% c(tmp_sitename,sites_to_pool[[i]])))

      new_ch_df <- rbind(new_ch_df_dropsites,removingsites_df) %>%
        dplyr::distinct(id,site,time,removed) %>%
        dplyr::group_by(id,site,time) %>%
        dplyr::summarize(id = dplyr::first(id),
                  site = dplyr::first(site),
                  time = dplyr::first(time),
                  removed = dplyr::last(removed))

      check <- new_ch_df %>%
        dplyr::group_by(id,site) %>%
        dplyr::mutate(dups = dplyr::n()) %>%
        dplyr::filter(dups > 1)

      if (nrow(check) > 0) message("Issue with pooling")

    }

  }



  observed_relative_min_max <-
    list(obs_min_time = obs_min_time,
         obs_max_time = obs_max_time,
         obs_min_a = obs_min_a,
         obs_max_a = obs_max_a,
         max_a_overall = max(set_max_a), # max_age_class
         set_max_a = set_max_a,
         set_min_a = set_min_a,
         # orig_max_a = max_a, # original max_a [but excluding pooled sites]
         max_rel_age_class = max(obs_max_a) - obs_min_a + 1,
         obsaux_min_a = obsaux_min_a,
         obsaux_max_a = obsaux_max_a,
         setaux_max_a = obsaux_max_a - min(obs_min_a) + 1
    )

  # if (any(obs_max_a > s4t_config$obs_max_a)) {
  #   message(paste0("Observed max ages exceed max ages in s4t_config:\n",
  #                  "Obs: ",paste0(obs_max_a,collapse = " "),
  #                  "\nSet: "),paste0(s4t_config$obs_max_a,collapse = " "))
  # }
  #
  # if (any(obs_min_a > s4t_config$obs_min_a)) {
  #   message(paste0("Observed min ages exceed max ages in s4t_config:\n",
  #                  "Obs: ",paste0(obs_min_a,collapse = " "),
  #                  "\nSet: "),paste0(s4t_config$obs_min_a,collapse = " "))
  # }






  sites_order <- s4t_config$sites_order


  tmp_removed_df <- new_ch_df %>%
    dplyr::ungroup() %>%
    dplyr::select(id,removed) %>%
    # mutate(removed = ifelse(removed == 1))
    dplyr::distinct(id,removed) %>%
    dplyr::mutate(site = removed) %>%
    dplyr::left_join(data.frame(site = s4t_config$sites_names,
                                order = 1:length(s4t_config$sites_names)),
                     by = "site") %>%
    dplyr::arrange(order) %>%
    dplyr::mutate(order = ifelse(is.na(order),0,order),
           removed = order) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(removed = dplyr::first(removed)) %>%
    dplyr::select(id,removed) # removed is now the number site


  wide_new_ch_df <- new_ch_df %>%
    dplyr::left_join(sites_order,by = "site") %>%
    dplyr::arrange(order) %>%
    dplyr::select(id,site,time) %>%
    dplyr::group_by(id) %>%
    tidyr::pivot_wider(names_from = site,
                values_from = time,
                values_fill = 0) %>%
    dplyr::left_join(tmp_removed_df, by = "id")

  ids_inboth <- intersect(aux_age_df_abbr$id,wide_new_ch_df$id)
  ids_not_in_ch <- setdiff(aux_age_df_abbr$id,wide_new_ch_df$id)
  ids_not_in_aux <- setdiff(wide_new_ch_df$id,aux_age_df_abbr$id)

  if (length(ids_not_in_ch) > 0) {
    message(paste0("Note, there are IDs in aux_age_df that are not in ch_df, n = ",length(ids_not_in_ch)))
  }

  if (length(ids_not_in_aux) > 0) {
    message(paste0("Dropping IDs in ch_df that are not in aux_age_df, n = ",length(ids_not_in_aux)))
  }


  keep_id_df <- data.frame(id = ids_inboth, keep_ids = TRUE)


  new_ch_aux_df <- wide_new_ch_df %>%
    dplyr::left_join(keep_id_df, by = "id") %>%
    dplyr::filter(!is.na(keep_ids)) %>%
    dplyr::select(-keep_ids) %>%
    dplyr::left_join(aux_age_df_abbr, by = "id")

  # new_ch_aux_df <- wide_new_ch_df %>%
  #   dplyr::filter(id %in% ids_inboth) %>%
  #   dplyr::left_join(aux_age_df_abbr, by = "id")

  # drop id column
  obs_ch <- as.matrix(new_ch_aux_df[,2:ncol(wide_new_ch_df)])
  obs_aux <- as.data.frame(new_ch_aux_df[,c((ncol(wide_new_ch_df) + 1):ncol(new_ch_aux_df))])


  # keep all auxiliary info
  all_aux <- aux_age_df_abbr %>%
    dplyr::select(-id)



  ### need to perform some checks first

  # is each individual recorded at each site once (not more than once)

  suppressMessages(repeatobservations <- new_ch_df %>%
                     dplyr::left_join(keep_id_df, by = "id") %>%
                     dplyr::filter(!is.na(keep_ids)) %>%
                     dplyr::select(-keep_ids) %>%
                     # dplyr::filter(id %in% ids_inboth) %>%
                     dplyr::group_by(id,site) %>%
                     dplyr::summarize(id_site = dplyr::n()) %>%
                     dplyr::filter(id_site > 1))

  suppressMessages(timedifferenceincaptures <- new_ch_df %>%
    dplyr::left_join(data.frame(site = s4t_config$sites_names,obs_max_a = obs_max_a,
                                obs_min_a = obs_min_a)) %>%
      dplyr::left_join(keep_id_df, by = "id") %>%
      dplyr::filter(!is.na(keep_ids)) %>%
      dplyr::select(-keep_ids) %>%
    # dplyr::filter(id %in% ids_inboth) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(first_obs = min_not_zero(time),
              last_obs = max(time),
              diff_time_obs = last_obs - first_obs) %>%
    dplyr::filter(diff_time_obs > obs_max_a - obs_min_a))


  suppressMessages(tmp_summary_dat <- new_ch_df %>%
                     dplyr::left_join(keep_id_df, by = "id") %>%
                     dplyr::filter(!is.na(keep_ids)) %>%
                     dplyr::select(-keep_ids) %>%
    # dplyr::filter(id %in% ids_inboth) %>%
    dplyr::group_by(site,time) %>%
    dplyr::summarize(observations = dplyr::n()) %>%
    dplyr::mutate(time = factor(time)))

  fewcaptures_insitetime <- tmp_summary_dat %>%
    dplyr::filter(observations < 10)

  suppressMessages(nocaptures_insitetime <- expand.grid(site = unique(new_ch_df$site),
                                       time = unique(new_ch_df$site)) %>%
    dplyr::left_join(tmp_summary_dat) %>%
    dplyr::filter(observations == 0))

  # checking that no removed individuals are encountered again

  # tmp_removed_ind <- obs_ch[obs_ch[,ncol(obs_ch)] != 0,] # 1:(ncol(obs_ch) - 1)
  # # tmp_removed_ind[2,3] <-7
  #
  # zombies <- (apply(tmp_removed_ind,MARGIN = 1,function(x) any(which(x[1:(length(x)-1)] != 0) > x[length(x)]))  )
  # zombie_obs <- tmp_removed_ind[zombies,]


  suppressWarnings(suppressMessages(zombie_individuals <- new_ch_df %>%
                                      dplyr::left_join(keep_id_df, by = "id") %>%
                                      dplyr::filter(!is.na(keep_ids)) %>%
                                      dplyr::select(-keep_ids) %>%
    # dplyr::filter(id %in% ids_inboth) %>%
    dplyr::filter(site %in% sites_order$site) %>%
    dplyr::left_join(sites_order, by = "site") %>%
    dplyr::group_by(id) %>%

    dplyr::filter(any(removed != FALSE)) %>%

    dplyr::arrange(id,order) %>%
    dplyr::mutate(removed_site = ifelse(removed != "FALSE",order,NA)) %>%
      # filter(id == "3D9.1C2D85B3A8")
    tidyr::fill(removed_site,.direction = "updown") %>%
    dplyr::filter((any(min(removed_site) < order,na.rm = TRUE)))
    ))

  ## time travelers

  # assume that site order is correct
  reversemovement <- new_ch_df %>%
    dplyr::left_join(sites_order, by = "site") %>%
    dplyr::arrange(id,order) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(!all(time == sort(time)))



  suppressMessages(max_obs_age_knownagefish <- new_ch_df %>%
    dplyr::left_join(aux_age_df_abbr[,c("id","ageclass","obs_time")], by = "id") %>%
    dplyr::group_by(id) %>%
    dplyr::filter(any(!is.na(ageclass))) %>%
    tidyr::fill(ageclass,obs_time) %>%
    dplyr::left_join(data.frame(site = s4t_config$sites_names,set_max_a = set_max_a)) %>%

    dplyr::mutate(first_obs = min_not_zero(time),
              last_obs = max(time),
              age_at_obs = dplyr::first(ageclass),
              obs_time = dplyr::first(obs_time),
              obs_time_max_a = obs_time + (min(obs_min_time) - 1) + set_max_a - age_at_obs,
              # set_obs_time = dplyr::first(obs_time)  - obs_min_time + 1,
              age_at_time = age_at_obs + (time - obs_time)#,
              # diff_time_obs = last_obs - (obs_time - obs_min_time + 1),
              # in_dataobs_max_age = obs_age + diff_time_obs
              ) %>%
    dplyr::filter(age_at_time > set_max_a)
  )

  init_relsite <- which(colSums(s4t_config$sites_config) == 0)
  init_relsite_names <- names(init_relsite)


  suppressMessages(missing_init_release_site <- new_ch_df %>%
                     dplyr::mutate(init_relsite = ifelse(site %in% init_relsite_names,TRUE,FALSE)) %>%
    dplyr::arrange(id,order) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(!any(init_relsite)))


  site_path <- s4t_config$site_path
  holdover_config <- s4t_config$holdover_config
  sites_config <- s4t_config$sites_config
  rowSums(s4t_config$holdover_config)

  # holdover_config[site_path[[1]],site_path[[1]]]

  # sapply(site_path,function(x) (as.matrix(holdover_config[x,x])))

  tmp_holdover_info <- sapply(site_path,function(x) colSums(as.matrix(holdover_config[x,x])))

  can_holdover <- matrix(FALSE,nrow = nrow(holdover_config),
                         ncol = ncol(holdover_config),
                         dimnames = list(rownames(holdover_config),rownames(holdover_config)))

  for (j in 1:(nrow(can_holdover) - 1)) {
    for (k in 2:ncol(can_holdover)) {
      if (k %in% site_path[[j]]) {
        k_index <- which(site_path[[j]] == k); k
        tmp <- sum(tmp_holdover_info[[j]][1:k_index])
        # print(paste0("j: ",j," k: ",k," holdover == ",tmp))
        can_holdover[j,k] <- tmp > 0
      }
    }
  }
  can_holdover_df <- can_holdover %>%
    as.data.frame() %>%
    dplyr::mutate(Site_j = rownames(can_holdover)) %>%
    tidyr::pivot_longer(cols = 1:ncol(can_holdover),names_to = "Site_k",values_to = "can_holdover")

  suppressMessages(identify_potential_holdovers <- new_ch_df %>% # by recap site
                     dplyr::left_join(data.frame(site = s4t_config$sites_names,holdover = colSums(s4t_config$holdover_config))) %>%
                     dplyr::left_join(sites_order, by = "site") %>%
                     dplyr::group_by(id) %>%
                     dplyr::arrange(id,order) %>%
                     dplyr::mutate(time_diff = time - dplyr::lag(time)) %>%
                     dplyr::mutate(time_diff = ifelse(is.na(time_diff),0,time_diff)) %>%
                     dplyr::filter(any(time_diff > 0 & holdover == 0))
                   )

  problematic_holdovers <- identify_potential_holdovers %>%
    dplyr::mutate(Site_j = dplyr::lag(site),
                  Site_k = site) %>%
    dplyr::filter(!is.na(Site_j)) %>%
    dplyr::left_join(can_holdover_df, by = c("Site_j","Site_k")) %>%
    dplyr::filter(can_holdover == FALSE) %>%
    dplyr::mutate(time = time + (obs_min_time - 1))



  potential_error_log <- list(repeatobservations = repeatobservations,
                              timedifferenceincaptures = timedifferenceincaptures,
                              fewcaptures_insitetime = fewcaptures_insitetime,
                              nocaptures_insitetime = nocaptures_insitetime,
                              zombie_individuals = zombie_individuals,
                              reversemovement = reversemovement,
                              max_obs_age_knownagefish = max_obs_age_knownagefish,
                              problematic_holdovers = problematic_holdovers,
                              missing_init_release_site = missing_init_release_site)

  message(paste0("\nError log:"))

  message(paste0("\nRepeat encounters at same site N = ",nrow(repeatobservations)))



  message(paste0("Individuals observed after being removed ('zombies') N = ",
                 length(unique(zombie_individuals$id))))

  message(paste0("Gap in observation times that exceed max difference in ages N = ",nrow(timedifferenceincaptures)))

  message(paste0("Holdovers observed between sites with only direct transitions N = ",nrow(problematic_holdovers)))

  # problematic_holdovers

  # reversemovement
  message(paste0("Reverse movements N = ",length(unique(reversemovement$id))))

  message(paste0("Known age individuals with ages greater than max age N = ",
                 length(unique(max_obs_age_knownagefish$id))))

  # missing_init_release_site

  message(paste0("Individuals with missing initial release site N = ",
                 length(unique(missing_init_release_site$id))))

  message("\nPotential errors:")

  message(paste0("Site/time combinations with 0 observations N = ",nrow(nocaptures_insitetime)))

  message(paste0("Site/time combinations with less than 10 observations N = ",nrow(fewcaptures_insitetime)))

  new_s4t_ch(obs_ch = obs_ch, # need to add something for removed individuals. maybe an additional column to obs_ch?
             obs_aux = obs_aux,
             all_aux = all_aux,
             ch_df = ch_df,
             set_max_a = set_max_a,
             s4t_config = s4t_config,
             observed_relative_min_max = observed_relative_min_max,
             potential_error_log = potential_error_log,
             cov_p = cov_p,
             cov_theta = cov_theta,
             call = call) # need to add stuff to relate age class and time to actual ages and years

  # potential_error_log
  # site name stuff

}

# to fix the "no visible binding for global variable" note:
site <- id <- time <- removed <- dups <- order <- id_site <-
  last_obs <- first_obs <- diff_time_obs <- observations <-
  removed_site <- ageclass <- obs_time <- obs_max_age <- time <-
  obs_age <- in_dataobs_max_age <- NULL


marginalize_ch <- function(s4t_ch) {
  # s4t_ch_all$ch$m_matrix
  # I think I need to drop all unused columns from auxilarly data to do this

  m_matrix_marg <- s4t_ch$ch$m_matrix %>%
    as.data.frame() %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    as.matrix()

  l_matrix_marg <- s4t_ch$ch$l_matrix %>%
    as.data.frame() %>%
    dplyr::group_by_all() %>%
    dplyr::summarize(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    as.matrix()




  # duplicated(as.data.frame(s4t_ch$ch$m_matrix))

  s4t_ch_marg <- s4t_ch

  s4t_ch_marg$ch$m_matrix_marg <- m_matrix_marg
  s4t_ch_marg$ch$m_matrix_larg <- l_matrix_marg

  class(s4t_ch_marg) <- "s4t_ch_marg"

  return(s4t_ch_marg)
}



#' Clean the capture history of a space4time capture history object
#'
#' @description
#' A short description...
#'
#' @param s4t_ch a `s4t_ch` object.
#' @returns a `clean_s4t_ch` object.
#' @details
#' Additional details...
#'
#'
#' @export
clean_s4t_ch_obs <- function(s4t_ch) {
  # s4t_ch$potential_error_log$timedifferenceincaptures

  # repeatobservations = repeatobservations,
  # timedifferenceincaptures = timedifferenceincaptures,
  # fewcaptures_insitetime = fewcaptures_insitetime,
  # nocaptures_insitetime = nocaptures_insitetime,
  # zombie_individuals = zombie_individuals,
  # reversemovement = reversemovement,
  # max_obs_age_knownagefish = max_obs_age_knownagefish

  ch_df <- as.data.frame(s4t_ch$obs_data$ch_df) %>%
    dplyr::mutate(drop_obs = FALSE)
  all_aux <- s4t_ch$obs_data$all_aux

  ## repeat observations

  message("Note: Cleaning method for repeat observations not implemented yet")

  identify_repeat_obs <- ch_df %>%
    dplyr::group_by(id,site) %>%
    dplyr::mutate(N = dplyr::n()) %>%
    dplyr::filter(N > 1) %>%
    dplyr::arrange(by = time) %>%
    dplyr::summarize(first_site_time_obs = dplyr::first(time))

  repeat_ids <- unique(identify_repeat_obs$id)
  if (length(repeat_ids) > 0 ) {
    # message("Note: Cleaning method for repeat observations is untested")

    for (i in 1:nrow(identify_repeat_obs)) {
      ch_df[ch_df$id == identify_repeat_obs$id[i] &
              ch_df$site == identify_repeat_obs$site[i] &
              ch_df$time > identify_repeat_obs$first_site_time_obs[i] ,]

      ch_df$drop_obs[ch_df$id == identify_repeat_obs$id[i] &
                       ch_df$site == identify_repeat_obs$site[i] &
                       ch_df$time > identify_repeat_obs$first_site_time_obs[i]] <- TRUE
    }
  }




  ## timedifferenceincaptures
  timedifferenceincaptures <- s4t_ch$potential_error_log$timedifferenceincaptures
  if (nrow(timedifferenceincaptures) > 0) {
    for (i in 1:nrow(timedifferenceincaptures)) {
      tmp_ind <- ch_df[ch_df$id == timedifferenceincaptures$id[i],]; tmp_ind
      first_obs <- min(tmp_ind$time); first_obs
      time_diff <- timedifferenceincaptures$diff_time_obs[i]; time_diff

#
      ch_df[ch_df$id == timedifferenceincaptures$id[i] &
              (ch_df$time - first_obs) ==  time_diff,]

      # if the observation exceeds the max time difference, drop that observation
      ch_df$drop_obs[ch_df$id == timedifferenceincaptures$id[i] &
                       (ch_df$time - first_obs) ==  time_diff] <- TRUE



    }

  }

  ## zombie_individuals
  zombie_individuals <- s4t_ch$potential_error_log$zombie_individuals

  if (nrow(zombie_individuals) > 0) {
    zombie_ids <- unique(zombie_individuals$id)
    for (i in 1:length(zombie_ids)) {
      tmp_ind <- ch_df[ch_df$id == zombie_ids[i],]; tmp_ind

      info_ind <- as.data.frame(zombie_individuals[zombie_individuals$id == zombie_ids[i],]); info_ind

      removed_site <- dplyr::first(info_ind$removed[which(info_ind$removed != "FALSE")])

      removed_site_ord <- info_ind$order[info_ind$site == removed_site]

      sites_past_removal <- info_ind$site[info_ind$order > removed_site_ord]

      for (j in sites_past_removal) {
        if (j %in% names(s4t_ch$s4t_config$sites_to_pool)) {
          j <- c(j,s4t_ch$s4t_config$sites_to_pool[[j]])
        }
        ch_df[ch_df$id == zombie_ids[i] &
                ch_df$site %in% j,]

        ch_df$drop_obs[ch_df$id == zombie_ids[i] &
                         ch_df$site %in% j] <- TRUE

      }


    }

  }


  ## reversemovement
  reversemovement <- s4t_ch$potential_error_log$reversemovement

  if (nrow(reversemovement) > 0) {
    revmove_ids <- unique(reversemovement$id)
    for (i in 1:length(revmove_ids)) {
      tmp_ind <- ch_df[ch_df$id == revmove_ids[i],]; tmp_ind

      info_ind <- as.data.frame(reversemovement[reversemovement$id == revmove_ids[i],]); info_ind

      # info_ind %>%
      #   dplyr::arrange(by = time) %>%
      #
      #
      # out_of_order_site <-

      # identify the last time the individual was observed, and drop that observation
      last_obs_site <- tmp_ind[which(tmp_ind$time == max(tmp_ind$time)),"site"]; last_obs_site


      ch_df[ch_df$id == revmove_ids[i] &
              ch_df$site %in% last_obs_site &
            ch_df$time == max(tmp_ind$time),]

      # drop observations
      ch_df$drop_obs[ch_df$id == revmove_ids[i] &
                       ch_df$site %in% last_obs_site &
                       ch_df$time == max(tmp_ind$time)] <- TRUE



    }

  }


  ## max_obs_age_knownagefish
  max_obs_age_knownagefish <- s4t_ch$potential_error_log$max_obs_age_knownagefish

  if (nrow(max_obs_age_knownagefish) > 0) {
    for (i in 1:nrow(max_obs_age_knownagefish)) {
      tmp_ind <- ch_df[ch_df$id == max_obs_age_knownagefish$id[i],]; tmp_ind

      first_obs <- min(tmp_ind$time) - s4t_ch$ch_info$observed_relative_min_max$obs_min_time + 1; first_obs
      obs_time <- max_obs_age_knownagefish$obs_time[i]; obs_time
      tmp_set_max_a <-max_obs_age_knownagefish$set_max_a[i]
      tmp_ageclass <-max_obs_age_knownagefish$ageclass[i]


      #
      ch_df[ch_df$id == max_obs_age_knownagefish$id[i] &
              tmp_set_max_a < tmp_ageclass + (ch_df$time - obs_time),]

      # if the observation exceeds the max obs time for the observed age, drop that observation
      ch_df$drop_obs[ch_df$id == max_obs_age_knownagefish$id[i] &
                       tmp_set_max_a < tmp_ageclass + (ch_df$time - obs_time)] <- TRUE

    }

  }


  ## problematic_holdovers
  problematic_holdovers <- s4t_ch$potential_error_log$problematic_holdovers

  if (nrow(problematic_holdovers) > 0) {
    for (i in 1:nrow(problematic_holdovers)) {
      tmp_ind <- ch_df[ch_df$id == problematic_holdovers$id[i],]; tmp_ind

      problematic_holdovers[i,]


      tmp_site <- problematic_holdovers$site[i]; tmp_site

      if (tmp_site %in% names(s4t_ch$s4t_config$sites_to_pool)) {
        tmp_site <- c(tmp_site,s4t_ch$s4t_config$sites_to_pool[[tmp_site]])
      }

      ch_df[ch_df$id == problematic_holdovers$id[i] &
              ch_df$site %in% tmp_site &
              ch_df$time >= problematic_holdovers$time[i],]

      # if the observation exceeds the max time difference, drop that observation
      ch_df$drop_obs[ch_df$id == problematic_holdovers$id[i] &
                       ch_df$site %in% tmp_site &
                       ch_df$time >= problematic_holdovers$time[i]] <- TRUE



    }

  }



  ## missing_init_release_site
  missing_init_release_site <- s4t_ch$potential_error_log$missing_init_release_site

  drop_ids_missing_inits <- unique(missing_init_release_site$id)
  ch_df$drop_obs[ch_df$id %in% drop_ids_missing_inits] <- TRUE

  # if (nrow(problematic_holdovers) > 0) {
  #   for (i in 1:nrow(problematic_holdovers)) {
  #     tmp_ind <- ch_df[ch_df$id == problematic_holdovers$id[i],]; tmp_ind
  #
  #     problematic_holdovers[i,]
  #
  #
  #     tmp_site <- problematic_holdovers$site[i]; tmp_site
  #
  #     if (tmp_site %in% names(s4t_ch$s4t_config$sites_to_pool)) {
  #       tmp_site <- c(tmp_site,s4t_ch$s4t_config$sites_to_pool[[tmp_site]])
  #     }
  #
  #     ch_df[ch_df$id == problematic_holdovers$id[i] &
  #             ch_df$site %in% tmp_site &
  #             ch_df$time >= problematic_holdovers$time[i],]
  #
  #     # if the observation exceeds the max time difference, drop that observation
  #     ch_df$drop_obs[ch_df$id == problematic_holdovers$id[i] &
  #                      ch_df$site %in% tmp_site &
  #                      ch_df$time >= problematic_holdovers$time[i]] <- TRUE
  #
  #
  #
  #   }
  #
  # }

  dropped_ch_df <- ch_df %>%
    as.data.frame() %>%
    dplyr::filter(drop_obs == TRUE) %>%
    dplyr::select(-drop_obs)

  cleaned_ch_df <- ch_df %>%
    as.data.frame() %>%
    dplyr::filter(drop_obs == FALSE) %>%
    dplyr::select(-drop_obs)

  dropped_ind <- length(setdiff(dropped_ch_df$id,cleaned_ch_df$id))

  message(paste0("N = ",nrow(dropped_ch_df)," observations and N = ",dropped_ind," individuals were dropped."))



  # if fewcaptures_insitetime or nocaptures_insitetime
  #  are greater than 1, send message that these cannot be resolved

  clean_s4t_ch <- list(cleaned_ch_df = cleaned_ch_df,dropped_ch_df = dropped_ch_df,
                       intermediate_ch_df = ch_df)
  class(clean_s4t_ch) <- "clean_s4t_ch"
  return(clean_s4t_ch)
}

n <- N <- drop_obs <- age_at_obs <- age_at_time <-
  holdover <- Site_j <- time_diff <- NULL

# validate_s4t_ch <- function(s4t_ch) {
#   # not up to date
#   if (sum(names(s4t_ch) != c("ch",
#                              "obs_data",
#                              "user_defined",
#                              "ch_info")) > 0) return(FALSE)
#
#   if (sum(names(s4t_ch$ch) != c("m_matrix",
#                                 "l_matrix",
#                                 "obs_lengthyear")) > 0) return(FALSE)
#
#   if (sum(names(s4t_ch$obs_data) != c("obs_ch",
#                                       "obs_lengthyear")) > 0) return(FALSE)
#
#   if (sum(names(s4t_ch$user_defined) != c("sites_config",
#                                           "holdover_config",
#                                           "max_a")) > 0) return(FALSE)
#
#   if (sum(names(s4t_ch$ch_info) != c('first_obs',
#                                      'n_batches',
#                                      'batches_list',
#                                      'first_sites',
#                                      'max_s_rel',
#                                      'max_t_recap',
#                                      'recap_sites',
#                                      'last_sites',
#                                      'recap_sites_not_last',
#                                      'observed_relative_min_max',
#                                      'potential_error_log')) > 0) return(FALSE)
#
#
#   return(TRUE)
# }
