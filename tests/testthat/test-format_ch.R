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

  site_arr <- linear_s4t_config(sites_names = 1:3,
                                holdover_sites = 1,
                                min_a =c(1,1,1),
                                max_a = c(3,3,3))

  expect_no_error( suppressMessages(s4t_ch(ch_df,
                               aux_age_df,
                               site_arr)
  ))

})


test_that("s4t_ch can be made from simulated data", {
  site_arr <- linear_s4t_config(sites_names = 1:3,
                                holdover_sites = 1,
                                min_a =c(1,1,1),
                                max_a = c(3,3,3))




  expect_no_error(suppressMessages(
    sim.dat <- sim_simple_s4t_ch(N = 2000,
                                 max_obs_year = 2
    )
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

  site_arr <- linear_s4t_config(sites_names = 1:3,
                                holdover_sites = 1,
                                min_a =c(1,1,1),
                                max_a = c(3,3,3))


 suppressMessages(ch <- s4t_ch(ch_df,
                   aux_age_df,
                   site_arr))
  expect_equal(nrow(ch$potential_error_log$reversemovement),3)


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

  site_arr <- linear_s4t_config(sites_names = 1:3,
                                holdover_sites = 1,
                                min_a =c(1,1,1),
                                max_a = c(3,3,3))


  suppressMessages(ch <- s4t_ch(ch_df,
                   aux_age_df,
                   site_arr)
  )
  expect_equal(nrow(ch$potential_error_log$timedifferenceincaptures),3)


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

  site_arr <- linear_s4t_config(sites_names = 1:3,
                                holdover_sites = 1,
                                min_a =c(1,1,1),
                                max_a = c(3,3,3))


  suppressMessages(ch <- s4t_ch(ch_df,
                   aux_age_df,
                   site_arr))
  expect_equal(nrow(ch$potential_error_log$max_obs_age_knownagefish),2)


})
