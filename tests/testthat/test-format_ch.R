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
