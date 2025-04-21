
# simulate_data_ppc <- function(s4t_ch,
#                               s4t_cjs_rstan,
#                               knownage = TRUE) {
#
#   if (knownage == FALSE) stop("Not yet implemented")
#
#   # maybe only use known age individuals for goodness of fit checks?
#
#   # obtain initial releases as a vector
#   releases = cbind(j = abc, # blah
#                    s = abc # blah
#   )
#
#   #obtain age length info stuff
#   age_info = cbind(obs_age, # = ... dropping unknown ages
#                    obs_time # = ... dropping unknown ages
#   )
#
#   group_info = abc # BLAH
#
#   # define true ch (i.e. true states)
#   sim_true_ch <- matrix(0,nrow = s4t_ch$abc,##,
#                         ncol = s4t_ch$abc)
#
#   # define observed ch (i.e. obs states)
#   sim_obs_ch <- matrix(0,nrow = s4t_ch$abc,##,
#                        ncol = s4t_ch$abc)
#
#   skip_i <- FALSE
#
#   for (i in 1:nrow(releases)) {
#     skip_i <- FALSE
#
#     first_j = releases[i,1]
#
#     for (j in site_path[[first_j]]) {
#       if (skip_i) next()
#
#       k = which(site_path[j,] == 1)
#       s = releases[i,2]
#       a1 = age_info[i,"obs_age"] + s - age_info[i,"obs_time"]
#       g = group_info[i]
#
#       # is sample faster?
#       a2 <- which(rmultinom(1,1,Theta[a1,,s,j,k,first_j,g])==1)
#
#       if (a2 == (max_a + 1)) {
#         skip_i <- TRUE
#         next() # this skips the rest of the capture history for this individual
#       }
#
#       t <- a2 - a1  + s
#
#       sim_true_ch[i,k] <- t
#
#       sim_obs_ch[i,k] <- ifelse(p_param[a1,a2,t,k,first_j,g],t,0) # HERE
#
#     }
#   }
#
#   # compute l-mat and m-mat
#   sim_m_matrix <- matrix(0,nrow = tot_entries_m, ncol = 8)
#   colnames(sim_m_matrix) <- c("j","k","s","t","b","g","obs_time","ageclass")
#   sim_l_matrix <- matrix(0,nrow = nrow(obs_ch),ncol = 6)
#   colnames(sim_l_matrix) <- c("j","s","b","g","obs_time","ageclass")
#
#
#   # placeholder, so 1 for now
#   sim_m_matrix[,"g"] <- 1 # blah
#   sim_l_matrix[,"g"] <- 1 # blah
#
#
#
#
#
#   counter <- 1
#   for (i in 1:nrow(sim_obs_ch)) {
#
#
#     obs <- which(sim_obs_ch[i,] != 0)
#     not_obs <- which(sim_obs_ch[i,] == 0)
#
#     first_rel <- obs[1]
#     last_rel <- obs[length(obs)]
#
#     #
#     entries_m <- length(obs) - 1
#     entries_l <- 1 # always 1
#
#     sim_l_matrix[i,1:2] <- c(last_rel,sim_obs_ch[i,last_rel])
#
#     # obtain location of first observation, which corresponds to the batch
#     sim_l_matrix[i,3] <- which(sim_obs_ch[i,]!= 0)[1]
#
#     sim_l_matrix[i,5] <- obs_aux[i,"obs_time"] # blah
#     sim_l_matrix[i,6] <- obs_aux[i,"ageclass"] # blah
#
#
#
#     if (entries_m > 0) {
#       for (z in 1:entries_m) {
#         sim_m_matrix[counter,1:4] <- c(obs[z],obs[z + 1],sim_obs_ch[i,obs[z]],sim_obs_ch[i,obs[z + 1]])
#         sim_m_matrix[counter,5] <- which(sim_obs_ch[i,]!= 0)[1]
#
#         sim_m_matrix[counter,7] <- obs_aux[i,"obs_time"] # blah
#         sim_m_matrix[counter,8] <- obs_aux[i,"ageclass"] # blah
#
#
#         counter <- counter + 1
#       }
#     }
#   }
#
#   # maybe actually format it like an m-array
#
#
#
#   # marginalize sim_m_matrix, sim_l_matrix, obs_m_matrix, and obs_l_matrix
#
#   # lambda_array[j,k,s,t,b,a1]
#   # chi_array[j,s,b,a1]
#
#   for (i in 1:nrow(sim_l_matrix)) {
#     j = sim_l_matrix[i,"j"]
#     s = sim_l_matrix[i,"s"]
#     b = sim_l_matrix[i,"b"]
#     g = sim_l_matrix[i,"g"]
#     a1 = sim_l_matrix[i,"a1"]
#
#     lambda_array[j,k,s,t,b,a1]
#
#
#
#   }
#
#   for (j in not_last_sites) {
#     for (a1 in min_a[j]:max_a[j]) {
#       for (s in min_rel_s[j]:max_rel_s[j]) {
#         for (b in 1:batchblah) {
#           for (g in 1:groupblah) {
#
#
#             sample_rows = (sim_obs_ch[,j] == s &
#                              a1 == (obs_age + s - obs_time) &
#                              b == obs_batch &
#                              g == group_infoblah)
#
#             sample_rows = (sim_obs_ch[,j] == s &
#                              a1 == (obs_age + s - obs_time)
#                            # b == obs_batch &
#                            # g == group_infoblah
#             )
#
#             num.released = sum(sample_rows)
#
#             sim_obs_ch[sample_rows,j:J]
#
#             next_obs_site <- j + apply(sim_obs_ch[sample_rows,(j+1):J],
#                                        MARGIN = 1, FUN = function(x) which(x != 0)[1])
#             next_obs_site <- ifelse(is.na(next_obs_site), 0,next_obs_site)
#
#             next_obs_time <- apply(sim_obs_ch[sample_rows,(j+1):J],
#                                    MARGIN = 1, FUN = function(x) x[which(x != 0)[1]])
#             next_obs_time <- ifelse(is.na(next_obs_time), 0,next_obs_time)
#
#             sim_vals = table(next_obs_site,next_obs_time) # add paired NAs with actual values
#             # to fill out the table
#
#             # also get obs_vals like sim_vals
#
#             max_t = s + max_a[next_site[k]] - a1
#             parameters = matrix(0, nrow = blah, ncol = max_a[next_site[k]])
#             parameters[1,1] <- chi_array[j,s,b,g,a1]
#             parameters[2:blah,2:blah] <- l_array[j,(j+1):J,s,s:max_t,b,a1]
#
#             exp_vals <- parameters * num.released
#
#             sim_FT <- sim_FT + sum((sqrt(sim_vals) - sqrt(exp_vals))^2)
#
#             obs_FT <- sim_FT + sum((sqrt(obs_vals) - sqrt(exp_vals))^2)
#
#           }
#         }
#       }
#     }
#     k = abc # blah
#     for (idk in idek) {
#       abc
#     }
#   }
#
#   lambda_array[]
#
#   # rowSums(Theta[1,1,,])
#   # Theta[1,1,,]
#
#
#
# }
