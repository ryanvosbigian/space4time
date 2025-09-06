
#' Simulate simple capture history
#'
#' @description
#' Simulates a small capture history with 3 sites, where there can be holdovers
#'     between the first 1st and 2nd sites.
#'
#' @param N An `integer` of the number of individuals
#' @param max_obs_year an `integer` of the number of total time intervals
#' @param prop_missing_age a `numeric` value between `0` and `1` of the proportion of
#'     individuals that do not have known ages. Default is 0.50.
#' @returns a `list` containing capture history `data.frame`, `list` containing the
#'     parameters.
#'
#' @export
sim_simple_s4t_ch <- function(N = 500,
                          max_obs_year = 2,
                          prop_missing_age = 0.5) {

  dat <- expand.grid(reps = 1:N,
                     obs_year = NA,
                     obs_age = NA,
                     obs_length = NA,
                     rst = 1,
                     lgr = NA,
                     bon = NA)


  dat$obs_length <- round(stats::runif(N,-6,6),1)#rnorm(nrow(dat))

  theta <- c(-Inf,-2,2,Inf)
  beta <- c(1)

  agelengthdat <- dat[,c("obs_age","obs_length")]

  agelengthdat$theta_0 <- stats::plogis(theta[1] - beta*agelengthdat$obs_length)
  agelengthdat$theta_1 <- stats::plogis(theta[2] - beta*agelengthdat$obs_length)
  agelengthdat$theta_2 <- stats::plogis(theta[3] - beta*agelengthdat$obs_length)
  agelengthdat$theta_3 <- stats::plogis(theta[4] - beta*agelengthdat$obs_length)

  agelengthdat$pi_3 <- agelengthdat$theta_3 - agelengthdat$theta_2
  agelengthdat$pi_2 <- agelengthdat$theta_2 - agelengthdat$theta_1
  agelengthdat$pi_1 <- agelengthdat$theta_1

  agelengthdat$pi_1 + agelengthdat$pi_2 + agelengthdat$pi_3

  dat$obs_age <- agelengthdat$age <- apply(agelengthdat[,c("pi_1","pi_2","pi_3")],MARGIN = 1,
                                           FUN = function(pi) sample(1:3,1,prob = pi))




  dat$obs_year <- sample(1:max_obs_year,nrow(dat),replace = TRUE)
  # dat$obs_age <- sample(1:3,nrow(dat),replace = TRUE)
  dat$broodyear <- dat$obs_year - dat$obs_age

  true_ch <- matrix(0,nrow = nrow(dat),ncol = 3)
  # true_ch[,1] <- sample(1:2,size = N,replace = TRUE)

  obs_ch <- matrix(0,nrow = nrow(dat),ncol = 3)
  obs_ch[,1] <- true_ch[,1] <- dat$obs_year
  broodyear <- dat$broodyear #+ 1 - min(dat$broodyear)

  BY_i <- dat$broodyear + 1 - min(dat$broodyear)

  obs_length <- dat$obs_length

  obs_lengthyear <- cbind(obs_year = dat$obs_year, obs_length = obs_length,obs_age = dat$obs_age)


  Theta <- array(data = 0, dim = c(s = max_obs_year +2,j = 2,a1=3,a2=4))
  for (y in 1:max_obs_year) {
    for (a1 in 1:3) {
      tmp_theta <- stats::runif(3,0.1,0.9)
      Theta[y,1,a1,a1] <- tmp_theta[a1]
      for (a2 in (a1+1):3) {
        if (a1 == 3) next()
        Theta[y,1,a1,a2] <- tmp_theta[a2] * prod(1 - tmp_theta[a1:(a2-1)])
      }
      Theta[y,1,a1,4] <- prod((1 - tmp_theta[a1:3]))

    }
  }

  diag(Theta[1,2,,]) <- stats::runif(3,min = 0.2,max = 0.8)
  diag(Theta[2,2,,]) <- stats::runif(3,min = 0.2,max = 0.8)
  diag(Theta[3,2,,]) <- stats::runif(3,min = 0.2,max = 0.8)
  diag(Theta[4,2,,]) <- stats::runif(3,min = 0.2,max = 0.8)

  Theta[1,2,,4] <- 1 - diag(Theta[1,2,,])
  Theta[2,2,,4] <- 1 - diag(Theta[2,2,,])
  Theta[3,2,,4] <- 1 - diag(Theta[3,2,,])
  Theta[4,2,,4] <- 1 - diag(Theta[4,2,,])


  p_prob <- matrix(stats::runif(4 * 2,0.5,0.95),
                   nrow = 4,
                   ncol = 2)

  rowSums(Theta[1,1,,])
  Theta[1,1,,]

  true_ch
  for (i in 1:N) {
    ## to second site
    a1 <- obs_lengthyear[i,"obs_age"]
    a2 <- which(stats::rmultinom(1,1,Theta[true_ch[i,1],1,a1,])==1)
    s <- true_ch[i,1]
    t <- a2 - a1  + true_ch[i,1]
    true_ch[i,2] <- ifelse(a2 ==4, 0,t) # if dead, record as dead, else record new age

    if (true_ch[i,2] == 0) next()

    obs_ch[i,2] <- ifelse(stats::rbinom(1,1,(true_ch[i,2] != 0) * p_prob[t,1]),
                          true_ch[i,2],
                          0)


    ## to third site
    a1 <- obs_lengthyear[i,"obs_age"] + true_ch[i,2] - obs_lengthyear[i,"obs_year"]
    a2 <- which(stats::rmultinom(1,1,Theta[true_ch[i,2],2,a1,])==1)
    s <- true_ch[i,2]
    t <- a2 - a1  + true_ch[i,2]
    true_ch[i,3] <- ifelse(a2 ==4, 0,t) # if dead, record as dead, else record new age

    if (true_ch[i,3] == 0) next()

    obs_ch[i,3] <- ifelse(stats::rbinom(1,1,(true_ch[i,3] != 0) * p_prob[t,2]),
                          true_ch[i,3],
                          0)

  }

  colnames(obs_lengthyear) <- c("OBS_YR","FL","OBS_AGE")

  obs_lengthyear[,"OBS_AGE"] <- ifelse(stats::rbinom(N,1,prop_missing_age),
                                       obs_lengthyear[,"OBS_AGE"],
                                       NA)


  overall_surv <- matrix(NA,nrow = 1 * max_obs_year * 3,ncol = 4)
  colnames(overall_surv) <- c("j","s","a1","survival")

  Theta
  counter <- 1
  for (j in 1) {
    for (s in 1:max_obs_year) {
      for (a1 in 1:3) {
        overall_surv[counter,] <- c(j,s,a1,sum(Theta[s,j,a1,1:3]))

        counter <- counter + 1
      }
    }
  }

  cohort_surv <- matrix(NA,nrow = 2 * (max_obs_year+2) * 3*3,ncol = 5)
  colnames(cohort_surv) <- c("j","s","a1","a2","survival")

  Theta
  counter <- 1
  for (j in 1:2) {
    for (s in 1:c(max_obs_year,max_obs_year + 2)[j]) {
      for (a1 in 1:3) {
        for (a2 in a1:3) {
          cohort_surv[counter,] <- c(j,s,a1,a2,Theta[s,j,a1,a2])

          counter <- counter + 1
        }

      }
    }
  }

  cohort_surv <- cohort_surv[!is.na(cohort_surv[,1]),]

  obs_ch <- cbind(id = 1:nrow(obs_ch), obs_ch)
  colnames(obs_ch) <- c("id", 1:3)

  #
  obs_lengthyear <- cbind(id = 1:nrow(obs_lengthyear), obs_lengthyear)
  colnames(obs_lengthyear) <- c("id", "obs_time", "FL", "ageclass")

  ch_df <- obs_ch %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = 2:4,
                        names_to = "site",
                        values_to = "time") %>%
    dplyr::mutate(removed = FALSE) %>%
    dplyr::filter(time != 0)


  s4t_config <- linear_s4t_config(sites_names = 1:3,
                                holdover_sites = 1,
                                min_a =c(1,1,1),
                                max_a = c(3,3,3))

  suppressMessages(s4t_ch <- s4t_ch(ch_df = ch_df,
                   aux_age_df = obs_lengthyear,
                   s4t_config = s4t_config
  ))

  return(list(s4t_ch = s4t_ch,
              params = list(Theta = Theta,
                            p_prob = p_prob,
                            overall_surv = overall_surv,
                            cohort_surv = cohort_surv)))



}


#' #' Simulate simple capture history with 4 sites
#' #'
#' #' @description
#' #' Simulates a small capture history with 4 sites, where there can be holdovers
#' #'     betweeen the first 1st and 2nd sites and the 2nd and 3rd sites.
#' #'
#' #' @param N An `integer` of the number of individuals
#' #' @param max_obs_year an `integer` of the number of total time intervals
#' #' @param prop_missing_age a `numeric` value between `0` and `1` of the proportion of
#' #'     individuals that do not have known ages.
#' #' @returns a `list` containing capture history `data.frame`, `list` containing the
#' #'     params.
#' #'
#' #'
#' sim_2_release_site_s4t_ch <- function(N = 500,
#'                                           max_obs_year = 2,
#'                                           prop_missing_age = 0.5) {
#'
#'   dat <- expand.grid(reps = 1:N,
#'                      obs_year = NA,
#'                      obs_age = NA,
#'                      obs_length = NA,
#'                      obs_site = NA,
#'                      rst = 1,
#'                      lgr = NA,
#'                      bon = NA)
#'
#'
#'   dat$obs_length <- round(stats::runif(N,-6,6),1)#rnorm(nrow(dat))
#'
#'   theta <- c(-Inf,-2,2,Inf)
#'   beta <- c(1)
#'
#'   agelengthdat <- dat[,c("obs_age","obs_length")]
#'
#'   agelengthdat$theta_0 <- stats::plogis(theta[1] - beta*agelengthdat$obs_length)
#'   agelengthdat$theta_1 <- stats::plogis(theta[2] - beta*agelengthdat$obs_length)
#'   agelengthdat$theta_2 <- stats::plogis(theta[3] - beta*agelengthdat$obs_length)
#'   agelengthdat$theta_3 <- stats::plogis(theta[4] - beta*agelengthdat$obs_length)
#'
#'   agelengthdat$pi_3 <- agelengthdat$theta_3 - agelengthdat$theta_2
#'   agelengthdat$pi_2 <- agelengthdat$theta_2 - agelengthdat$theta_1
#'   agelengthdat$pi_1 <- agelengthdat$theta_1
#'
#'   agelengthdat$pi_1 + agelengthdat$pi_2 + agelengthdat$pi_3
#'
#'   dat$obs_age <- agelengthdat$age <- apply(agelengthdat[,c("pi_1","pi_2","pi_3")],MARGIN = 1,
#'                                            FUN = function(pi) sample(1:3,1,prob = pi))
#'
#'
#'
#'
#'   dat$obs_year <- sample(1:max_obs_year,nrow(dat),replace = TRUE)
#'   dat$obs_site <- sample(1:2,nrow(dat),replace = TRUE)
#'   # dat$obs_age <- sample(1:3,nrow(dat),replace = TRUE)
#'   dat$broodyear <- dat$obs_year - dat$obs_age
#'
#'   true_ch <- matrix(0,nrow = nrow(dat),ncol = 4)
#'   # true_ch[,1] <- sample(1:2,size = N,replace = TRUE)
#'
#'   obs_ch <- matrix(0,nrow = nrow(dat),ncol = 4)
#'   for (i in 1:N) {
#'     obs_ch[i,dat$obs_site[i]] <- true_ch[i,dat$obs_site[i]] <- dat$obs_year[i]
#'   }
#'
#'   broodyear <- dat$broodyear #+ 1 - min(dat$broodyear)
#'
#'   BY_i <- dat$broodyear + 1 - min(dat$broodyear)
#'
#'   obs_length <- dat$obs_length
#'
#'   obs_lengthyear <- cbind(obs_year = dat$obs_year, obs_length = obs_length,obs_age = dat$obs_age,
#'                           obs_site = dat$obs_site)
#'
#'
#'   Theta <- array(data = 0, dim = c(s = max_obs_year +4,j = 3,a1=3,a2=4))
#'   for (j in 1:2) {
#'     for (y in 1:c(max_obs_year,max_obs_year + 2)[j]) {
#'       for (a1 in 1:3) {
#'         tmp_theta <- stats::runif(3,0.1,0.9)
#'         Theta[y,j,a1,a1] <- tmp_theta[a1]
#'         for (a2 in (a1+1):3) {
#'           if (a1 == 3) next()
#'           Theta[y,j,a1,a2] <- tmp_theta[a2] * prod(1 - tmp_theta[a1:(a2-1)])
#'         }
#'         Theta[y,j,a1,4] <- prod((1 - tmp_theta[a1:3]))
#'
#'       }
#'     }
#'   }
#'
#'
#'   diag(Theta[1,3,,]) <- stats::runif(3,min = 0.2,max = 0.8)
#'   diag(Theta[2,3,,]) <- stats::runif(3,min = 0.2,max = 0.8)
#'   diag(Theta[3,3,,]) <- stats::runif(3,min = 0.2,max = 0.8)
#'   diag(Theta[4,3,,]) <- stats::runif(3,min = 0.2,max = 0.8)
#'
#'   Theta[1,3,,4] <- 1 - diag(Theta[1,3,,])
#'   Theta[2,3,,4] <- 1 - diag(Theta[2,3,,])
#'   Theta[3,3,,4] <- 1 - diag(Theta[3,3,,])
#'   Theta[4,3,,4] <- 1 - diag(Theta[4,3,,])
#'
#'
#'   p_prob <- matrix(stats::runif(4 * 3,0.5,0.95),
#'                    nrow = 4,
#'                    ncol = 3)
#'
#'   # rowSums(Theta[1,1,,])
#'   # Theta[1,1,,]
#'
#'   true_ch
#'   for (i in 1:N) {
#'     ## to second site
#'     if (true_ch[i,1] != 0) {
#'       a1 <- obs_lengthyear[i,"obs_age"]
#'       a2 <- which(stats::rmultinom(1,1,Theta[true_ch[i,1],1,a1,])==1)
#'       s <- true_ch[i,1]
#'       t <- a2 - a1  + true_ch[i,1]
#'       true_ch[i,2] <- ifelse(a2 ==4, 0,t) # if dead, record as dead, else record new age
#'
#'       if (true_ch[i,2] == 0) next()
#'
#'       obs_ch[i,2] <- ifelse(stats::rbinom(1,1,(true_ch[i,2] != 0) * p_prob[t,1]),
#'                             true_ch[i,2],
#'                             0)
#'     }
#'
#'
#'
#'     ## to third site
#'     a1 <- obs_lengthyear[i,"obs_age"] + true_ch[i,2] - obs_lengthyear[i,"obs_year"]
#'     a2 <- which(stats::rmultinom(1,1,Theta[true_ch[i,2],2,a1,])==1)
#'     s <- true_ch[i,2]
#'     t <- a2 - a1  + true_ch[i,2]
#'     true_ch[i,3] <- ifelse(a2 ==4, 0,t) # if dead, record as dead, else record new age
#'
#'     if (true_ch[i,3] == 0) next()
#'
#'     obs_ch[i,3] <- ifelse(stats::rbinom(1,1,(true_ch[i,3] != 0) * p_prob[t,2]),
#'                           true_ch[i,3],
#'                           0)
#'
#'
#'     ## to third site
#'     a1 <- obs_lengthyear[i,"obs_age"] + true_ch[i,3] - obs_lengthyear[i,"obs_year"]
#'     a2 <- which(stats::rmultinom(1,1,Theta[true_ch[i,3],3,a1,])==1)
#'     s <- true_ch[i,3]
#'     t <- a2 - a1  + true_ch[i,3]
#'     true_ch[i,4] <- ifelse(a2 ==4, 0,t) # if dead, record as dead, else record new age
#'
#'     if (true_ch[i,4] == 0) next()
#'
#'     obs_ch[i,4] <- ifelse(stats::rbinom(1,1,(true_ch[i,4] != 0) * p_prob[t,2]),
#'                           true_ch[i,4],
#'                           0)
#'
#'   }
#'
#'   colnames(obs_lengthyear) <- c("OBS_YR","FL","OBS_AGE","OBS_SITE")
#'
#'   obs_lengthyear[,"OBS_AGE"] <- ifelse(stats::rbinom(N,1,1-prop_missing_age),
#'                                        obs_lengthyear[,"OBS_AGE"],
#'                                        NA)
#'
#'
#'
#'   overall_surv <- matrix(NA,nrow = 2 * (max_obs_year+2) * 3,ncol = 4)
#'   colnames(overall_surv) <- c("j","s","a1","survival")
#'
#'   Theta
#'   counter <- 1
#'   for (j in 1:2) {
#'     for (s in 1:c(max_obs_year,max_obs_year + 2)[j]) {
#'       for (a1 in 1:3) {
#'         overall_surv[counter,] <- c(j,s,a1,sum(Theta[s,j,a1,1:3]))
#'
#'         counter <- counter + 1
#'       }
#'     }
#'   }
#'
#'   overall_surv <- overall_surv[!is.na(overall_surv[,1]),]
#'
#'
#'   cohort_surv <- matrix(NA,nrow = 2 * (max_obs_year+2) * 3*3,ncol = 5)
#'   colnames(cohort_surv) <- c("j","s","a1","a2","survival")
#'
#'   Theta
#'   counter <- 1
#'   for (j in 1:2) {
#'     for (s in 1:c(max_obs_year,max_obs_year + 2)[j]) {
#'       for (a1 in 1:3) {
#'         for (a2 in a1:3) {
#'           cohort_surv[counter,] <- c(j,s,a1,a2,Theta[s,j,a1,a2])
#'
#'           counter <- counter + 1
#'         }
#'
#'       }
#'     }
#'   }
#'
#'   cohort_surv <- cohort_surv[!is.na(cohort_surv[,1]),]
#'
#' #
#' #   sites_config <- matrix(c(0, 1, 0, 0,
#' #                            0, 0, 1, 0,
#' #                            0, 0, 0, 1,
#' #                            0, 0, 0, 0), nrow = 4, byrow = TRUE)
#' #   colnames(sites_config) <- rownames(sites_config) <- 1:4
#' #
#' #   holdover_config <- matrix(c(0, 1, 0, 0,
#' #                               0, 0, 1, 0,
#' #                               0, 0, 0, 0,
#' #                               0, 0, 0, 0), nrow = 4, byrow = TRUE)
#' #
#' #   colnames(holdover_config) <- rownames(holdover_config) <- 1:4
#'
#'   s4t_config <- linear_s4t_config(sites_names = 1:4,
#'                                   holdover_sites = c(1,2),
#'                                   min_a =c(1,1,1,1),
#'                                   max_a = c(3,3,3,3))
#'
#'
#'   ch_df <- obs_ch
#'
#'   suppressMessages(s4t_ch <- s4t_ch(ch_df = ch_df,
#'     aux_age_df = obs_lengthyear,
#'     s4t_config = s4t_config
#'   ))
#'
#'   return(list(s4t_ch = s4t_ch,
#'               params = list(Theta = Theta,
#'                             p_prob = p_prob,
#'                             overall_surv = overall_surv,
#'                             cohort_surv = cohort_surv)))
#'
#'   # return(list(obs_ch = obs_ch,true_ch = true_ch,
#'   #             aux_age_df = broodyear,
#'   #             obs_lengthyear = obs_lengthyear,
#'   #             s4t_config = s4t_config,
#'   #             params = list(Theta = Theta,
#'   #                           p_prob = p_prob,
#'   #                           overall_surv = overall_surv,
#'   #                           cohort_surv = cohort_surv)
#'   # )
#'   # )
#'
#'
#' }
#'
