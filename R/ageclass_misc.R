


ageclass_call <- function(age_formula = ~ 1,
                          obs_aux,
                          # max_a, # drop this?
                          ll = TRUE) { # drop ll = TRUE? doesn't do anything?
  obs_aux <- as.data.frame(obs_aux)



  if ("ageclass" %in% colnames(obs_aux)) {
    obsageclass <- obs_aux[,"ageclass"]
  } else {
    obsageclass <- NULL
    if (ll == TRUE) {
      stop("Missing ageclass column in obs_aux")
    }
  }

  obs_aux_df <- as.data.frame(obs_aux)

  # obs_aux_df$obs_time <- as.factor(obs_aux_df$obs_time)

  mf_reg <- stats::model.frame(age_formula,data = as.data.frame(obs_aux))
  mt_reg <- attr(mf_reg,"terms")

  mod_mat_a_beta <- stats::model.matrix(mt_reg, mf_reg)

  tmp <- paste0("a_beta_",colnames(mod_mat_a_beta))
  tmp2 <- rep(0.1,length(tmp)); names(tmp2) <- tmp

  par_beta <- c(tmp2[-1])

  #
  ## Alpha

  new_max_a <- max(obs_aux[,"ageclass"],na.rm = TRUE)

  par_alpha <- rep(1,new_max_a - 1)
  par_alpha[1] <- -1
  names(par_alpha) <- paste0("a_alpha_",1:(new_max_a - 1))


  age_par <- c(par_alpha,par_beta)

  list(age_par = age_par,
       mod_mat_a_beta = mod_mat_a_beta,
       obsageclass = obsageclass)

}

ageclass_nll <- function(par,
                         max_a,
                         mod_mat_a_beta,
                         obsageclass = NULL,
                         ll = TRUE) {
  # made a change to to max_a so that it will be more coherent


  # low level function
  if (is.null(obsageclass) | mean(is.na(obsageclass)) == 1) {
    # new_max_a <- max_a
    new_max_a <- 1 + sum(grepl("a_alpha_",names(par)))
  } else {
    new_max_a <- max(obsageclass,na.rm = TRUE)
  }

  alpha <- vector()
  beta <- vector()



  for (a in 1:(new_max_a-1)) {
    alpha[a] <- par[paste0("a_alpha_",a)]
  }
  # This is probably bad practice. Could exponentiate it
  alpha[2:(new_max_a-1)] <- abs(alpha[2:(new_max_a-1)])

  # beta <- par[grepl("^a_beta",par)]

  # all beta parameters, just don't start with a_alpha ??
  beta <- c(0,par[grepl("^a_beta",names(par))])

  # beta[1] <- par["beta"]

  eta <- vector()
  eta[1] <- -Inf

  for (i in 1:(new_max_a - 1)) {
    eta[i+1] <- sum(alpha[1:i])
  }

  eta[new_max_a+1] <- Inf


  theta <- matrix(NA,nrow = nrow(mod_mat_a_beta),ncol = new_max_a+1)
  pi <- matrix(NA,nrow = nrow(mod_mat_a_beta),ncol = max_a)


  theta[,1] <- 0
  theta[,new_max_a + 1] <- 1
  for (i in 2:(new_max_a)) {
    mu <- (eta[i] - mod_mat_a_beta %*% beta)

    theta[,i] <- stats::plogis(mu)
  }

  pi[,1] <-  theta[,2]
  for (i in 2:new_max_a) {
    pi[,i] <- theta[,i+1] -  theta[,i]
  }

  if (new_max_a < max_a) {
    for (j in (new_max_a + 1):max_a) {
      pi[,j] <- 0
    }
  }

  if (sum(pi < 0) > 0) {
    utils::head(pi); utils::tail(pi)
    print(par)
    stop("negative pi")
  }



  if (ll) {

    alk_ll <- -sum(sapply(1:nrow(mod_mat_a_beta),
                          FUN = function(i) log(pi[i,obsageclass[i]])))

    if (is.na(alk_ll) | is.infinite(alk_ll)) {
      print(par)
      stop("Negative or NA nLL in agemodel")
    }

    return(alk_ll)

  } else{

    return(pi)
  }

}



alk_likelihood <- function(par,
                           obs_aux,
                           max_a,
                           ll = TRUE) {
  # simple age-length key. Multinomial logistic link ordinal regression.
  # Only one parameter for relating length to age.
  # Should make one where the effect of length varies by broodyear
  obs_aux <- as.data.frame(obs_aux)

  alpha <- vector()
  beta <- vector()

  for (a in 1:(max_a-1)) {
    alpha[a] <- par[paste0("alpha[",a,"]")]
  }

  beta[1] <- par["beta"]

  eta <- vector()
  eta[1] <- -Inf

  for (i in 1:(max_a - 1)) {
    eta[i+1] <- sum(alpha[1:i])
  }

  eta[max_a+1] <- Inf


  theta <- matrix(NA,nrow = nrow(obs_aux),ncol = max_a+1)
  pi <- matrix(NA,nrow = nrow(obs_aux),ncol = max_a)

  for (i in 1:(max_a + 1)) {
    theta[,i] <- stats::plogis(eta[i] - beta[1]*obs_aux[,"FL"] )
  }

  pi[,1] <-  theta[,2]
  for (i in 2:max_a) {
    pi[,i] <- theta[,i+1] -  theta[,i]
  }


  if (ll) {
    alk_ll <- -sum(sapply(1:nrow(obs_aux),FUN = function(i)
      stats::dmultinom(diag(max_a)[obs_aux[i,"ageclass"],],size = 1,
                prob = pi[i,],log = TRUE)))

    return(alk_ll)

  } else{

    return(pi)
  }
}
