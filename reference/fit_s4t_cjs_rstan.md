# Fit space-for-time mark-recapture model in Bayesian framework

Fits the time- and age-specific space-for-time mark-recapture model in a
Bayesian framework using Stan.

## Usage

``` r
fit_s4t_cjs_rstan(
  p_formula,
  theta_formula,
  ageclass_formula,
  groups = NULL,
  s4t_ch,
  chains = 3,
  warmup = 500,
  iter = 1000,
  fixed_age = TRUE,
  ...
)
```

## Arguments

- p_formula:

  a formula for detection probabilities using default indices (a1, a2,
  s, t, j, k, r, and g) or covariates.

- theta_formula:

  a formula for transition probabilities using default indices (a1, a2,
  s, t, j, k, r, and g) or covariates.

- ageclass_formula:

  a formula for the effect structure of the age-class sub-model.

- groups:

  a `character` vector containing the names of the covariates that
  comprise the groups.

- s4t_ch:

  a `s4t_ch` object

- chains:

  an `integer` of the number of chains to run.

- warmup:

  an `integer` of the number of warmup iterations.

- iter:

  an `integer` of the number of total (warmup + actual) iterations.

- fixed_age:

  a `logical` object that determines whether the ageclass model will be
  run as a separate model (TRUE) or whether it is estimated along with
  the CJS model (FALSE).

- ...:

  further arguents to pass to
  [`rstan::sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

## Value

a `s4t_cjs_rstan` object.

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_s4t_cjs_rstan(p_formula = ~ t,
                        theta_formula = ~ a1 * a2 * s * j,
                        ageclass_formula = ageclass ~ FL,
                        fixed_age = TRUE,
                        s4t_ch = sim.dat$s4t_ch)
} # }
```
