
test_that("models can fit ml", {
  skip_on_cran()

  set.seed(1)
  sim.dat <- sim_simple_s4t_ch(N = 2000)

  # suppressMessages(
  #   ch <- s4t_ch(
  #     ch_df = sim.dat$ch_df,
  #     aux_age_df = sim.dat$aux_age_df,
  #     s4t_config = sim.dat$s4t_config
  #   )
  # )

  expect_no_error(suppressMessages(
    m1.fixed <- fit_s4t_cjs_ml(
      p_formula = ~ t,
      theta_formula = ~ a1 * a2 * s * j,
      ageclass_formula = ~ FL,
      fixed_age = TRUE,
      s4t_ch = sim.dat$s4t_ch
    )
  ))

  expect_no_error(suppressMessages(
    m1.notfixed <- fit_s4t_cjs_ml(
      p_formula = ~ t,
      theta_formula = ~ a1 * a2 * s * j,
      ageclass_formula = ~ FL,
      fixed_age = FALSE,
      s4t_ch = sim.dat$s4t_ch
    )
  ))

})



test_that("models can fit rstan", {
  skip_on_cran()

  set.seed(1)
  sim.dat <- sim_simple_s4t_ch(N = 800)

  expect_no_error(suppressMessages(
    m1.fixed <- fit_s4t_cjs_rstan(
      p_formula = ~ t,
      theta_formula = ~ a1 * a2 * s * j,
      ageclass_formula = ~ FL,
      fixed_age = TRUE,
      s4t_ch = sim.dat$s4t_ch,
      chains = 2,
      warmup = 200,
      iter = 400
    )
  ))

  expect_no_error(suppressMessages(
    m1.notfixed <- fit_s4t_cjs_rstan(
      p_formula = ~ t,
      theta_formula = ~ a1 * a2 * s * j,
      ageclass_formula = ~ FL,
      fixed_age = FALSE,
      s4t_ch = sim.dat$s4t_ch,
      chains = 2,
      warmup = 200,
      iter = 400
    )
  ))

})
