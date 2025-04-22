test_that("s4t_ch is made", {

  ch_df <- data.frame(id = c(1,1,1,
                             2,2,
                             3,3,
                             4,
                             5,
                             6),
                      site = c(1,2,3,
                               1,2,
                               1,3,
                               1,
                               1,
                               1),
                      time = c(1,3,3,
                               2,3,
                               1,3,
                               2,
                               1,
                               1),
                      removed = c(FALSE,FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,
                                  FALSE,
                                  FALSE
                                  ))

  aux_age_df <- data.frame(id = 1:6,
                           obs_site = rep(1,6),
                           ageclass = c(1,2,1,1,2,1),
                           obs_time = c(1,2,1,2,1,1),
                           Covariate1 = c(3,1,2,1,2,1))


  expect_no_error(  s4t_cjs_ch(ch_df,
                               aux_age_df,
                               min_a = c(1,1,1),
                               max_a = c(3,3,3),
                               sites_names = 1:3,
                               sites_config = matrix(c(0,1,0,
                                                       0,0,1,
                                                       0,0,0),
                                                     nrow = 3,
                                                     ncol = 3,
                                                     byrow = TRUE,
                                                     dimnames = list(1:3,1:3)),
                               holdover_config = matrix(c(0,1,0,
                                                          0,0,0,
                                                          0,0,0),
                                                        nrow = 3,
                                                        ncol = 3,
                                                        byrow = TRUE,
                                                        dimnames =  list(1:3,1:3))
  ))

})


test_that("s4t_ch can be made from simulated data", {
  sites_config <- matrix(c(0, 1, 0,
                           0, 0, 1,
                           0, 0, 0), nrow = 3, byrow = TRUE)
  colnames(sites_config) <- rownames(sites_config) <- 1:3

  holdover_config <- matrix(c(0, 1, 0,
                              0, 0, 0,
                              0, 0, 0), nrow = 3, byrow = TRUE)

  colnames(holdover_config) <- rownames(holdover_config) <- 1:3

  sim.dat <- simulate_data(N = 1000,
                           max_obs_year = 2
  )
  # sim.dat$obs_ch <- cbind(id = 1:nrow(sim.dat$obs_ch), sim.dat$obs_ch)
  # colnames(sim.dat$obs_ch) <- c("id", 1:3)
  #
  # sim.dat$obs_lengthyear <- cbind(id = 1:nrow(sim.dat$obs_lengthyear), sim.dat$obs_lengthyear)
  # colnames(sim.dat$obs_lengthyear) <- c("id", "obs_time", "FL", "ageclass")


  # ch_df <- sim.dat$obs_ch %>%
  #   as.data.frame() %>%
  #   tidyr::pivot_longer(cols = 2:4,
  #                names_to = "site",
  #                values_to = "time") %>%
  #   dplyr::mutate(removed = FALSE) %>%
  #   dplyr::filter(time != 0)

  expect_no_error(s4t_cjs_ch(
    ch_df = sim.dat$ch_df,
    aux_age_df = sim.dat$obs_lengthyear,
    max_a = c(3,3,3),
    min_a = c(1,1,1),
    sites_names = as.character(1:3),
    sites_config = sites_config,
    holdover_config = holdover_config
  ))

})


test_that("error is thrown by bad capture history - reverse movement", {

  ch_df <- data.frame(id = c(1,1,1,
                             2,2,
                             3,3,
                             4,
                             5,
                             6),
                      site = c(1,2,3,
                               1,2,
                               1,3,
                               1,
                               1,
                               1),
                      time = c(1,3,1,
                               2,3,
                               1,3,
                               2,
                               1,
                               1),
                      removed = c(FALSE,FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,
                                  FALSE,
                                  FALSE
                      ))

  aux_age_df <- data.frame(id = 1:6,
                           obs_site = rep(1,6),
                           ageclass = c(1,2,1,1,2,1),
                           obs_time = c(1,2,1,2,1,1),
                           Covariate1 = c(3,1,2,1,2,1))

  ch <- s4t_cjs_ch(ch_df,
                   aux_age_df,
                   min_a = c(1,1,1),
                   max_a = c(3,3,3),
                   sites_names = 1:3,
                   sites_config = matrix(c(0,1,0,
                                           0,0,1,
                                           0,0,0),
                                         nrow = 3,
                                         ncol = 3,
                                         byrow = TRUE,
                                         dimnames = list(1:3,1:3)),
                   holdover_config = matrix(c(0,1,0,
                                              0,0,0,
                                              0,0,0),
                                            nrow = 3,
                                            ncol = 3,
                                            byrow = TRUE,
                                            dimnames =  list(1:3,1:3)))
  expect_equal(nrow(ch$ch_info$potential_error_log$reversemovement),3)


})


test_that("error is thrown by bad capture history - exceed gap in obs time", {

  ch_df <- data.frame(id = c(1,1,1,
                             2,2,
                             3,3,
                             4,
                             5,
                             6),
                      site = c(1,2,3,
                               1,2,
                               1,3,
                               1,
                               1,
                               1),
                      time = c(1,5,5,
                               2,3,
                               1,3,
                               2,
                               1,
                               1),
                      removed = c(FALSE,FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,
                                  FALSE,
                                  FALSE
                      ))

  aux_age_df <- data.frame(id = 1:6,
                           obs_site = rep(1,6),
                           ageclass = c(NA,2,1,1,2,1),
                           obs_time = c(1,2,1,2,1,1),
                           Covariate1 = c(3,1,2,1,2,1))

  ch <- s4t_cjs_ch(ch_df,
                   aux_age_df,
                   min_a = c(1,1,1),
                   max_a = c(3,3,3),
                   sites_names = 1:3,
                   sites_config = matrix(c(0,1,0,
                                           0,0,1,
                                           0,0,0),
                                         nrow = 3,
                                         ncol = 3,
                                         byrow = TRUE,
                                         dimnames = list(1:3,1:3)),
                   holdover_config = matrix(c(0,1,0,
                                              0,0,0,
                                              0,0,0),
                                            nrow = 3,
                                            ncol = 3,
                                            byrow = TRUE,
                                            dimnames =  list(1:3,1:3)))
  expect_equal(nrow(ch$ch_info$potential_error_log$timedifferenceincaptures),3)


})


test_that("error is thrown by bad capture history - exceed max age", {

  ch_df <- data.frame(id = c(1,1,1,
                             2,2,
                             3,3,
                             4,
                             5,
                             6),
                      site = c(1,2,3,
                               1,2,
                               1,3,
                               1,
                               1,
                               1),
                      time = c(1,3,3,
                               2,3,
                               1,3,
                               2,
                               1,
                               1),
                      removed = c(FALSE,FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,FALSE,
                                  FALSE,
                                  FALSE,
                                  FALSE
                      ))

  aux_age_df <- data.frame(id = 1:6,
                           obs_site = rep(1,6),
                           ageclass = c(3,2,1,1,2,1),
                           obs_time = c(1,2,1,2,1,1),
                           Covariate1 = c(3,1,2,1,2,1))

  ch <- s4t_cjs_ch(ch_df,
                   aux_age_df,
                   min_a = c(1,1,1),
                   max_a = c(3,3,3),
                   sites_names = 1:3,
                   sites_config = matrix(c(0,1,0,
                                           0,0,1,
                                           0,0,0),
                                         nrow = 3,
                                         ncol = 3,
                                         byrow = TRUE,
                                         dimnames = list(1:3,1:3)),
                   holdover_config = matrix(c(0,1,0,
                                              0,0,0,
                                              0,0,0),
                                            nrow = 3,
                                            ncol = 3,
                                            byrow = TRUE,
                                            dimnames =  list(1:3,1:3)))
  expect_equal(nrow(ch$ch_info$potential_error_log$max_obs_age_knownagefish),3)


})
