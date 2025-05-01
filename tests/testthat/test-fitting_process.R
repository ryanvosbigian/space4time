
test_that("models can fit ml ", {
  skip_on_cran()

  set.seed(1)
  sim.dat <- simulate_data(N = 2000)

  suppressMessages(
    ch <- s4t_cjs_ch(
      ch_df = sim.dat$ch_df,
      aux_age_df = sim.dat$aux_age_df,
      min_a = sim.dat$min_a,
      max_a = sim.dat$max_a,
      sites_names = sim.dat$sites_names,
      sites_config = sim.dat$sites_config,
      holdover_config = sim.dat$holdover_config
    )
  )

  expect_no_error(suppressMessages(
    m1 <- fit_s4t_cjs_ml(
      p_formula = ~ t,
      theta_formula = ~ a1 * a2 * s * j,
      ageclass_formula = ~ FL,
      fixed_age = TRUE,
      s4t_ch = ch
    )
  ))


})



test_that("models can fit rstan", {
  skip_on_cran()

  set.seed(1)
  sim.dat <- simulate_data(N = 800)

  suppressMessages(
    ch <- s4t_cjs_ch(
      ch_df = sim.dat$ch_df,
      aux_age_df = sim.dat$aux_age_df,
      min_a = sim.dat$min_a,
      max_a = sim.dat$max_a,
      sites_names = sim.dat$sites_names,
      sites_config = sim.dat$sites_config,
      holdover_config = sim.dat$holdover_config
    )
  )

  expect_no_error(suppressMessages(
    m1.s <- fit_s4t_cjs_rstan(
      p_formula = ~ t,
      theta_formula = ~ a1 * a2 * s * j,
      ageclass_formula = ~ FL,
      fixed_age = TRUE,
      s4t_ch = ch,
      chains = 2,
      warmup = 200,
      iter = 400
    )
  ))

})
