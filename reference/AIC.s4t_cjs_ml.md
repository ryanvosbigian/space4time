# Compute AIC of fitted model objects

Compute AIC of one or more fitted `s4t_cjs_ml` model objects

## Usage

``` r
# S3 method for class 's4t_cjs_ml'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  A `s4t_cjs_ml` model object

- ...:

  Optionally more fitted model objects

- k, :

  The penalty parameter, taken to be 2. Not used but needed for generic
  consistency

## Value

If one object is provided, just returns a numeric value with the
corresponding AIC. If more than one is provided, it returns a
`data.frame` with rows corresponding to the objects and columns
representing the number of parameters estimated (`df`), and the AIC

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_s4t_cjs_ml(p_formula = ~ t,
                     theta_formula = ~ a1 * a2 * s * j,
                     ageclass_formula = ~ FL,
                     fixed_age = TRUE,
                     s4t_ch = sim.dat$s4t_ch)
m2 <- fit_s4t_cjs_ml(p_formula = ~ t * a1 * a2,
                     theta_formula = ~ a1 * a2 * s * j,
                     ageclass_formula = ~ FL,
                     fixed_age = TRUE,
                     s4t_ch = sim.dat$s4t_ch)
AIC(m1, m2)
} # }
```
