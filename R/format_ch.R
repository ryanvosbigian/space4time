
#' @importFrom magrittr %>%
process_ch <- function(obs_ch,obs_aux,removed_sites,sites = NULL,batches_list = NULL) {
  obs_ch <- as.matrix(obs_ch)
  # obs_aux <- as.matrix(obs_aux)
  # obs_ch <- sim.dat$obs_ch
  # obs_lengthyear <- sim.dat$obs_lengthyear
  max_t <- max(obs_ch)

  # (max_t + 1)*(max_t + 1 - 1) * (max_t + 1 - 2)

  tot_entries_m <- sum(apply(obs_ch, MARGIN = 1,FUN = function(x) sum(x != 0) - 1))

  additional_aux_var <- setdiff(colnames(obs_aux),c("id","obs_time","ageclass"))

  m_matrix <- matrix(0,nrow = tot_entries_m, ncol = 8)
  colnames(m_matrix) <- c("j","k","s","t","b","g","obs_time","ageclass")
  l_matrix <- matrix(0,nrow = nrow(obs_ch),ncol = 6)
  colnames(l_matrix) <- c("j","s","b","g","obs_time","ageclass")


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

    # obtain location of first observation, which corresponds to the batch
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


new_s4t_cjs_ch <- function(obs_ch,
                           obs_aux,
                           set_max_a,
                           sites_config,
                           holdover_config,
                           observed_relative_min_max,
                           potential_error_log,
                           call) {

  mat_obs_ch <- as.matrix(obs_ch[,colnames(sites_config)])
  removed_sites <- (obs_ch[,ncol(obs_ch)])
  # obs_aux <- as.matrix(obs_aux)



  ch <- process_ch(mat_obs_ch,obs_aux,removed_sites)
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
      # b_diff_batches <- setdiff(f,batches)
      # intersect(b_diff_batches,first_sites)

      if (length(f_int_1st) == length(f)){

        keep_going <- FALSE
      } else {
        tmp_f1 <- intersect(f,first_sites)
        tmp_f2 <- f# setdiff(b,tmp_b1)
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



  # first_obs <- apply(obs_ch,MARGIN = 1,FUN = function(x) which(x != 0)[1])
  first_obs <- apply(obs_ch[,1:(ncol(obs_ch)-1)],MARGIN = 1,FUN = function(x) which.min.not.zero(x)[1])

  # sort(unique(first_obs))

  # which(sim.dat$obs_ch[,1])

  n_batches <- length(sort(unique(first_obs)))# sum(colSums(sites) == 0)
  batches <- sort(unique(first_obs)) # which(colSums(sites) == 0)

  batches_list <- list()

  first_sites <- which(colSums(sites_config) == 0)

  for (j in 1:nrow(sites_config)) {

    b <- j
    # b <- c()

    keep_going <- TRUE
    while (keep_going) {
      b_int_batches <- intersect(b,batches)
      b_diff_batches <- setdiff(b,batches)

      # NEED to check whether all sites are properly explored

      explored <- TRUE
      for (i in 1:length(b)) {
        # if there are any sites in the path that haven't been explored.
        if (length(setdiff(site_1stin_path[[b[i]]],b)) > 0) {
          explored <- FALSE
        }
      }


      if (length(intersect(b,batches)) == length(b) &
          explored == TRUE){

        keep_going <- FALSE
      } else {
        # save all batches encountered
        tmp_b1 <- intersect(b,batches)
        tmp_b2 <- b# setdiff(b,tmp_b1)

        # save all the b's that are a first site.
        tmp_b3 <- setdiff(b,first_sites)
        tmp_b4 <- intersect(b,first_sites)

        # loop through all the b's that are not a first site to explore it further
        app <- c()
        if (length(tmp_b3) > 0) {
          for (i in 1:length(tmp_b3)) {
            app <-  c(app,which(sites_config[,tmp_b3[i]] == 1))
          }
        }

        b <- (sort(c(tmp_b1,tmp_b4,app))) # unique


      }


    }


    batches_list[[j]] <- b



  } # end for loop
  batches_list


  n_stations <- ncol(sites_config)

  # should have as many entries as n_stations
  max_s_rel <- rep(0,n_stations); names(max_s_rel) <- 1:n_stations
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

  # max_t for sites with recaps. Typically second to last n_stations (when there is only 1 release site)
  max_t_recap <- rep(0,n_stations); names(max_t_recap) <- 1:n_stations
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


  ## Determines the minimum and maximum age for each site-time.
  # This is needed, because otherwise, the likelihood will include
  # transitions between ages that are not possible.

  # determine first captures
  first_cap_j  <- apply(obs_ch[,1:(ncol(obs_ch)-1)],MARGIN = 1,
                        FUN = function(x) which(x != 0)[1])


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

        min_ageclass_mat[k,s] <- min(c(min_age_sk,
                                       existing_val_min_tmp),na.rm = TRUE)


        # max values are simple, just given
        max_ageclass_mat[k,s] <- set_max_a[k]


      }

    }
  }




  recap_sites <- which(colSums(sites_config) > 0)
  last_sites <- unique(unlist(lapply(site_path,FUN = max)))

  recap_sites_not_last <- recap_sites[!(recap_sites %in% last_sites)]


  not_last_sites <- c(1:n_stations)[-last_sites]


  s4t_ch <- list(ch = list(m_matrix = m_matrix,
                           l_matrix = l_matrix,
                           m_aux_df = m_aux_df,
                           l_aux_df = l_aux_df,
                           obs_aux = obs_aux),
                 obs_data = list(obs_ch = obs_ch,
                                 obs_aux = obs_aux),
                 call = call,
                 user_defined = list(sites_config = sites_config,
                                     holdover_config = holdover_config,
                                     max_a = observed_relative_min_max$orig_max_a,
                                     set_max_a = set_max_a,
                                     sites_names = rownames(sites_config)),
                 ch_info = list(first_obs = first_obs,
                                n_batches = n_batches,
                                batches_list = batches_list,
                                n_stations = n_stations,
                                first_sites = first_sites,
                                site_path = site_path,
                                max_s_rel = max_s_rel,
                                max_t_recap = max_t_recap,
                                set_max_a = set_max_a,
                                set_min_a = observed_relative_min_max$set_min_a,
                                recap_sites = recap_sites,
                                last_sites = last_sites,
                                recap_sites_not_last = recap_sites_not_last,
                                observed_relative_min_max = observed_relative_min_max,
                                potential_error_log = potential_error_log,
                                min_ageclass_mat = min_ageclass_mat,
                                max_ageclass_mat = max_ageclass_mat,
                                init_rel_j = init_rel_j,
                                init_rel_times = init_rel_times)
  )
  class(s4t_ch) = "s4t_cjs_ch"

  return(s4t_ch)
}


#' Create capture history object
#'
#' @description
#' A short description...
#'
#' @param ch_df `data.frame` object containing each capture. See LINK VIGNETTE and details.
#' @param aux_age_df `data.frame` containing auxiliary data for each individual. See LINK VIGNETTE and details.
#' @param cov_df `NULL` or `data.frame` of covariates aligned with `j,k,s,t,b,r` indices. See LINK VIGNETTE and details.
#' @param min_a a `vector` of the minimum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param max_a a `vector` of the maximum ageclass individuals in a site can be.
#'     Must be the same length and order as `sites_names`
#' @param sites_names a character `vector` of the site names
#' @param sites_config a `matrix` that determine how sites are linked together.
#'     Must be in the same order as sites_names. See LINK VIGNETTE and details.
#' @param holdover_config a `matrix` that determine whether individuals can holdover in between
#'      sites. In the same format as sites_configs. See LINK VIGNETTE and details.
#' @param sites_to_pool a list of character vectors that contain the names of sites
#'      to be pooled together and treated as one site. See LINK VIGNETTE and details.
#' @details
#' Additional details...
#' @examples
#'
#' ch_df <- data.frame(id = c(1,1,1,
#'                            2,2,
#'                            3,3,
#'                            4,
#'                            5,
#'                            6),
#'                      site = c(1,2,3,
#'                               1,2,
#'                               1,3,
#'                               1,
#'                               1,
#'                               1),
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
#'                           obs_site = rep(1,6),
#'                           ageclass = c(1,2,1,1,2,1),
#'                           obs_time = c(1,2,1,2,1,1),
#'                           Covariate1 = c(3,1,2,1,2,1))
#'
#' ch <- s4t_cjs_ch(ch_df,aux_age_df,
#'                  min_a = c(1,1,1),
#'                  max_a = c(3,3,3),
#'                  sites_names = 1:3,
#'                  sites_config = matrix(c(0,1,0,
#'                                          0,0,1,
#'                                          0,0,0),
#'                                          nrow = 3,
#'                                          ncol = 3,
#'                                          byrow = TRUE,
#'                                          dimnames = list(1:3,1:3)),
#'                 holdover_config = matrix(c(0,1,0,
#'                                            0,0,0,
#'                                            0,0,0),
#'                                            nrow = 3,
#'                                            ncol = 3,
#'                                            byrow = TRUE,
#'                                            dimnames = list(1:3,1:3)),
#'                 )
#'
#'
#'
#'
#' @export
s4t_cjs_ch <- function(ch_df,
                       aux_age_df,
                       cov_df = NULL,
                       min_a,
                       max_a,
                       sites_names,
                       sites_config,
                       holdover_config,
                       sites_to_pool = NULL) {

  call <- list(ch_df = ch_df,
               aux_age_df = aux_age_df,
               cov_df = cov_df,
               min_a = min_a,
               max_a = max_a,
               sites_names = sites_names,
               sites_config = sites_config,
               holdover_config = holdover_config,
               sites_to_pool = sites_to_pool)

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


  if (!is.null(cov_df)) {
    if (!("data.frame" %in% class(cov_df))) {
      stop("cov_df must be a data.frame")
    }

  }


  max_a <- as.integer(max_a)
  min_a <- as.integer(min_a)

  if (sum(is.na(max_a))>1) {
    stop("max_a must be coercible to an integer")
  }

  if (sum(is.na(min_a))>1) {
    stop("min_a must be coercible to an integer")
  }

  if (length(max_a) != ncol(sites_config)) {
    stop(paste0("max_a should be length = ",ncol(sites_config)))
  }

  if (length(min_a) != ncol(sites_config)) {
    stop(paste0("minx_a should be length = ",ncol(sites_config)))
  }


  sites_names <- as.character(sites_names)
  if (sum(is.na(sites_names))>1) {
    stop("sites_names must be coercible to a character vector")
  }


  if (!is(sites_config,"matrix") | length(dim(sites_config)) != 2) {
    stop("sites_config must be a matrix")
  }


  if (!is(holdover_config,"matrix") | length(dim(holdover_config)) != 2) {
    stop("holdover_config must be a matrix")
  }


  if (!is.null(sites_to_pool)) {
    if (!is(sites_to_pool,"list")) {
      stop("sites_to_pool must be a list")
      # should also check if the names match
    }
    if (length(setdiff(names(sites_to_pool),sites_names)) > 0) {
      stop("Names of sites_to_pool not in sites_names")
    }
    for (i in 1:length(sites_to_pool)) {
      if (length(setdiff(sites_to_pool[[i]],sites_names)) > 0) {
        stop("Names of sites in sites_to_pool not in sites_names")
      }
    }
  }





  # check that sites_names and sites_config, holdover_config are the right dimensions
  if (length(unique(c(nrow(sites_config),ncol(sites_config),
                      nrow(holdover_config ),ncol(holdover_config ),
                      length(sites_names), length(max_a),
                      length(min_a)))) != 1) {
    stop("Dimensions of sites_config, holdover_config")
  }

  # check that if configs are named, that the names match sites_names
  if (!is.null(rownames(sites_config))) {
    if (any(!(rownames(sites_config) == sites_names))) {
      stop("If present, row and column names of sites_config must equal and in the same order as sites_names")
    }
  }
  if (!is.null(colnames(sites_config))) {
    if (any(!(colnames(sites_config) == sites_names))) {
      stop("If present, row and column names of sites_config must equal and in the same order as sites_names")
    }
  }
  if (!is.null(rownames(holdover_config))) {
    if (any(!(rownames(holdover_config) == sites_names))) {
      stop("If present, row and column names of holdover_config must equal and in the same order as sites_names")
    }
  }
  if (!is.null(colnames(holdover_config))) {
    if (any(!(colnames(holdover_config) == sites_names))) {
      stop("If present, row and column names of holdover_config must equal and in the same order as sites_names")
    }
  }



  # ch_df$id
  # ch_df$site

  # change type of ch_df$id and ch_df$site to characters
  ch_df$id <- as.character(ch_df$id)
  ch_df$site <- as.character(ch_df$site)






  tmp_sites_to_remove <- setdiff(unique(ch_df$site),sites_names)

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
  min_obs_time <- min_not_zero(c(ch_df$time,aux_age_df$obs_time))
  max_obs_time <- max(c(ch_df$time,aux_age_df$obs_time),na.rm = TRUE)



  # the minimum obs age-class will be age-class 1 (even if it is 0, 2, etc.)
  min_obs_age <- min(aux_age_df$ageclass,na.rm = TRUE)

  # min_obs_age <- min_not_zero(aux_age_df$ageclass)
  max_obs_age <- max(aux_age_df$ageclass,na.rm = TRUE)

  # sets the minimum age_class to 1.
  set_max_a <- max_a - min_obs_age + 1
  set_min_a <- min_a  - min_obs_age + 1

  max_rel_age_class = max(max_a) - min_obs_age + 1



  ##
  # age - auxiliary variables

  # aux_age_df
  aux_age_df_abbr <- aux_age_df
  aux_age_df_abbr$ageclass <- aux_age_df_abbr$ageclass - min_obs_age + 1
  aux_age_df_abbr$obs_time <- aux_age_df_abbr$obs_time - min_obs_time + 1


  ## next create obs_ch.

  new_ch_df <- data.frame(id = ch_df$id,
                          site = ch_df$site,
                          time = ch_df$time - min_obs_time + 1,
                          removed = as.logical(ch_df$removed))
  new_ch_df$removed <- as.character(ifelse(new_ch_df$removed,new_ch_df$site,FALSE))


  # pool sites

  using_sites <- sites_names

  if (!is.null(sites_to_pool)) {
    for (i in 1:length(sites_to_pool)) {

      tmp_sitename <- names(sites_to_pool)[i]

      if (tmp_sitename %in% sites_to_pool[[i]]) sites_to_pool[[i]] <- setdiff(sites_to_pool[[i]],tmp_sitename)


      ## ADD THE CHECK RIGHT HERE
      tmp_sconfig <- sites_config[which(sites_names == sites_to_pool[[i]] | sites_names == tmp_sitename),]
      tmp_seqal <- apply(tmp_sconfig,MARGIN = 1,FUN = function(x) length(unique(x)))

      tmp_hconfig <- holdover_config[which(sites_names == sites_to_pool[[i]] | sites_names == tmp_sitename),]
      tmp_heqal <- apply(tmp_hconfig,MARGIN = 1,FUN = function(x) length(unique(x)))

      if (any(tmp_seqal != 1) | any(tmp_heqal != 1)) {
        stop(paste0("Cannot pool sites with different site arrangements or holdovers: ",tmp_sitename," with: ",
                    paste0(sites_to_pool[[i]],collapse = ", ")))
      }

      # remove the sites that we are pooling
      using_sites <- setdiff(using_sites,sites_to_pool[[i]])

      using_sites <- ifelse(using_sites == tmp_sitename,paste0(tmp_sitename,"_plus"),using_sites)
      sites_names <- using_sites

      # sites_config and holdovers config to drop:
      keepthesesites <- !(sites_names %in% sites_to_pool[[i]])
      sites_config <- sites_config[keepthesesites,keepthesesites]
      holdover_config <- holdover_config[keepthesesites,keepthesesites]

      # dropping min_a and max_a values
      min_a <- min_a[keepthesesites]
      max_a <- max_a[keepthesesites]
      set_min_a <- set_min_a[keepthesesites]
      set_max_a <- set_max_a[keepthesesites]

      colnames(holdover_config) <- rownames(holdover_config) <- sites_names
      colnames(sites_config) <- rownames(sites_config) <- sites_names

      sites <- sites_config
      holdover_config <- holdover_config

      removingsites_df <- new_ch_df %>%
        dplyr::filter(site %in% sites_to_pool[[i]]) %>%
        dplyr::mutate(site = paste0(tmp_sitename,"_plus")) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(id) %>%
        dplyr::summarize(id = dplyr::first(id),
                  site = paste0(tmp_sitename,"_plus"),
                  time = min(time),
                  removed = (ifelse(any(removed != FALSE),TRUE,FALSE)))
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
        dplyr::summarize(dups = dplyr::n()) %>%
        dplyr::filter(dups > 1)

      if (nrow(check) > 0) message("Issue with pooling")

    }

  }


  observed_relative_min_max <-
    list(min_obs_time = min_obs_time,
         max_obs_time = max_obs_time,
         min_obs_age = min_obs_age,
         max_obs_age = max_obs_age,
         max_age_class = max(set_max_a),
         set_max_a = set_max_a,
         set_min_a = set_min_a,
         orig_max_a = max_a, # original max_a [but excluding pooled sites]
         max_rel_age_class = max(max_a) - min_obs_age + 1
    )









  site_order <- data.frame(site = colnames(sites_config),
                           ord = 1:ncol(sites_config))


  tmp_removed_df <- new_ch_df %>%
    dplyr::ungroup() %>%
    dplyr::select(id,removed) %>%
    # mutate(removed = ifelse(removed == 1))
    dplyr::distinct(id,removed) %>%
    dplyr::mutate(site = removed) %>%
    dplyr::left_join(site_order, by = "site") %>%
    dplyr::arrange(ord) %>%
    dplyr::mutate(ord = ifelse(is.na(ord),0,ord),
           removed = ord) %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(removed = dplyr::first(removed)) %>%
    dplyr::select(id,removed) # removed is now the number site


  wide_new_ch_df <- new_ch_df %>%
    dplyr::left_join(site_order,by = "site") %>%
    dplyr::arrange(ord) %>%
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
    message(paste0("Dropping IDs in aux_age_df that are not in ch_df, n = ",length(ids_not_in_ch)))
  }

  if (length(ids_not_in_aux) > 0) {
    message(paste0("Dropping IDs in ch_df that are not in aux_age_df, n = ",length(ids_not_in_aux)))
  }

  new_ch_aux_df <- wide_new_ch_df %>%
    dplyr::filter(id %in% ids_inboth) %>%
    dplyr::left_join(aux_age_df_abbr, by = "id")

  # drop id column
  obs_ch <- as.matrix(new_ch_aux_df[,2:ncol(wide_new_ch_df)])
  obs_aux <- as.data.frame(new_ch_aux_df[,c((ncol(wide_new_ch_df) + 1):ncol(new_ch_aux_df))])




  ### need to perform some checks first

  # is each individual recorded at each site once (not more than once)

  suppressMessages(repeatobservations <- ch_df %>%
                     dplyr::filter(id %in% ids_inboth) %>%
                     dplyr::group_by(id,site) %>%
                     dplyr::summarize(id_site = dplyr::n()) %>%
                     dplyr::filter(id_site > 1))

  suppressMessages(timedifferenceincaptures <- ch_df %>%
    dplyr::left_join(data.frame(site = sites_names,max_a = max_a)) %>%
    dplyr::filter(id %in% ids_inboth) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(first_obs = min_not_zero(time),
              last_obs = max(time),
              diff_time_obs = last_obs - first_obs) %>%
    dplyr::filter(diff_time_obs > max_a - min_obs_age))


  suppressMessages(tmp_summary_dat <- ch_df %>%
    dplyr::filter(id %in% ids_inboth) %>%
    dplyr::group_by(site,time) %>%
    dplyr::summarize(observations = dplyr::n()) %>%
    dplyr::mutate(time = factor(time)))

  fewcaptures_insitetime <- tmp_summary_dat %>%
    dplyr::filter(observations < 10)

  suppressMessages(nocaptures_insitetime <- expand.grid(site = unique(ch_df$site),
                                       time = unique(ch_df$site)) %>%
    dplyr::left_join(tmp_summary_dat) %>%
    dplyr::filter(observations == 0))

  # checking that no removed individuals are encountered again

  # tmp_removed_ind <- obs_ch[obs_ch[,ncol(obs_ch)] != 0,] # 1:(ncol(obs_ch) - 1)
  # # tmp_removed_ind[2,3] <-7
  #
  # zombies <- (apply(tmp_removed_ind,MARGIN = 1,function(x) any(which(x[1:(length(x)-1)] != 0) > x[length(x)]))  )
  # zombie_obs <- tmp_removed_ind[zombies,]


  suppressWarnings(suppressMessages(zombie_individuals <- ch_df %>%
    dplyr::filter(id %in% ids_inboth) %>%
    dplyr::filter(site %in% site_order$site) %>%
    dplyr::left_join(site_order, by = "site") %>%
    dplyr::group_by(id) %>%
    dplyr::filter(any(removed)) %>%
    dplyr::arrange(id,ord) %>%
    dplyr::mutate(removed_site = ifelse(removed,ord,NA)) %>%
    tidyr::fill(removed_site,.direction = "updown") %>%
    dplyr::filter((any(min(removed_site) < ord,na.rm = TRUE)))))

  ## time travelers

  # assume that site order is correct
  reversemovement <- ch_df %>%
    dplyr::left_join(site_order, by = "site") %>%
    dplyr::arrange(id,ord) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(!all(time == sort(time)))


  ## NEED TO REDO
  suppressMessages(max_obs_age_knownagefish <- ch_df %>%
    dplyr::left_join(aux_age_df, by = "id") %>%
    dplyr::left_join(data.frame(site = sites_names,max_a = max_a)) %>%
    dplyr::group_by(id) %>%
    dplyr::filter(any(!is.na(ageclass))) %>%
    tidyr::fill(ageclass,obs_time) %>%
    dplyr::mutate(first_obs = min_not_zero(time),
              last_obs = max(time),
              obs_age = dplyr::first(ageclass),
              obs_time = dplyr::first(obs_time),
              diff_time_obs = last_obs - obs_time,
              obs_max_age = obs_age + diff_time_obs) %>%
    dplyr::filter(obs_max_age > max_a))




  potential_error_log <- list(repeatobservations = repeatobservations,
                              timedifferenceincaptures = timedifferenceincaptures,
                              fewcaptures_insitetime = fewcaptures_insitetime,
                              nocaptures_insitetime = nocaptures_insitetime,
                              zombie_individuals = zombie_individuals,
                              reversemovement = reversemovement,
                              max_obs_age_knownagefish = max_obs_age_knownagefish)

  message(paste0("Potential error log:"))

  message(paste0("Number of individuals encountered more than once at a site: ",nrow(repeatobservations)))

  message(paste0("Number of individuals with a gap in observation times that exceed difference in minimum and maximum ages: ",nrow(timedifferenceincaptures)))

  message(paste0("Number of site/time combinations with less than 10 observations: ",nrow(fewcaptures_insitetime)))

  message(paste0("Number of site/time combinations with no observations: ",nrow(nocaptures_insitetime)))

  message(paste0("Number of 'zombies' (individuals observed after being removed): ",
                 length(unique(zombie_individuals$id))))

  # reversemovement
  message(paste0("Number of individuals with reverse movements: ",length(unique(reversemovement$id))))

  message(paste0("Number of known age individuals with observed ages greater than max_a: ",
                 length(unique(max_obs_age_knownagefish$id))))

  new_s4t_cjs_ch(obs_ch = obs_ch, # need to add something for removed individuals. maybe an additional column to obs_ch?
                 obs_aux = obs_aux,
                 set_max_a = set_max_a,
                 sites_config = sites_config,
                 holdover_config = holdover_config,
                 observed_relative_min_max = observed_relative_min_max,
                 potential_error_log = potential_error_log,
                 call = call) # need to add stuff to relate age class and time to actual ages and years

  # potential_error_log
  # site name stuff

}

# to fix the "no visible binding for global variable" note:
site <- id <- time <- removed <- dups <- ord <- id_site <-
  last_obs <- first_obs <- diff_time_obs <- observations <-
  removed_site <- ageclass <- obs_time <- obs_max_age <- time <-
  obs_age <- NULL


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

validate_s4t_ch <- function(s4t_ch) {
  # not up to date
  if (sum(names(s4t_ch) != c("ch",
                             "obs_data",
                             "user_defined",
                             "ch_info")) > 0) return(FALSE)

  if (sum(names(s4t_ch$ch) != c("m_matrix",
                                "l_matrix",
                                "obs_lengthyear")) > 0) return(FALSE)

  if (sum(names(s4t_ch$obs_data) != c("obs_ch",
                                      "obs_lengthyear")) > 0) return(FALSE)

  if (sum(names(s4t_ch$user_defined) != c("sites_config",
                                          "holdover_config",
                                          "max_a")) > 0) return(FALSE)

  if (sum(names(s4t_ch$ch_info) != c('first_obs',
                                     'n_batches',
                                     'batches_list',
                                     'first_sites',
                                     'max_s_rel',
                                     'max_t_recap',
                                     'recap_sites',
                                     'last_sites',
                                     'recap_sites_not_last',
                                     'observed_relative_min_max',
                                     'potential_error_log')) > 0) return(FALSE)


  return(TRUE)
}
