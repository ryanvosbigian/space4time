# Fit space-for-time mark-recapture model in likelihood framework

Fits the time- and age-specific space-for-time mark-recapture model
using maximum likelihood.

## Usage

``` r
fit_s4t_cjs_ml(
  p_formula,
  theta_formula,
  ageclass_formula,
  groups = NULL,
  s4t_ch,
  ndeps = 0.001,
  lmm = 5,
  maxit = 500,
  fixed_age = FALSE
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
  comprise the groups. Default is no groups (`groups = NULL`).

- s4t_ch:

  a `s4t_ch` object

- ndeps:

  a `numeric` value passed to `optim`. See `optim` documentation.

- lmm:

  an `integer` value passed to `optim`. See `optim` documentation.

- maxit:

  an `integer` value passed to optim. See `optim` documentation.

- fixed_age:

  a `logical` value that determines whether the age-class model will be
  run as a separate model (`TRUE`) or whether it is estimated along with
  the mark-recapture model (`FALSE`).

## Value

a `s4t_cjs_ml` object.

## Details

Additional details...

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_s4t_cjs_ml(p_formula = ~ t,
                     theta_formula = ~ a1 * a2 * s * j,
                     ageclass_formula = ~ FL,
                     fixed_age = TRUE,
                     s4t_ch = sim.dat$s4t_ch)
} # }


```
