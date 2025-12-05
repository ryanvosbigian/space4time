# Calculate variance-covariance matrix for a fitted model object

Calculate variance-covariance matrix for a `s4t_cjs_ml` model object.

## Usage

``` r
# S3 method for class 's4t_cjs_ml'
vcov(object, ...)
```

## Arguments

- object:

  a `s4t_cjs_ml` object fitted using
  [`fit_s4t_cjs_ml()`](https://ryanvosbigian.github.io/space4time/reference/fit_s4t_cjs_ml.md).

- ...:

  Other arguments. Not used (needed for generic consistency)

## Value

The variance-covariance matrix of coefficients from a fitted
`s4t_cjs_ml` object

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_s4t_cjs_ml(p_formula = ~ t,
                     theta_formula = ~ a1 * a2 * s * j,
                     ageclass_formula = ~ FL,
                     fixed_age = TRUE,
                     s4t_ch = sim.dat$s4t_ch)

vcov(m1)
} # }
```
