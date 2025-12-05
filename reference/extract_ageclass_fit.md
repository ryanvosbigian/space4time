# Extract `s4t_ageclass_model` object from a `s4t_cjs_rstan` or `s4t_cjs_ml` object

Extract `s4t_ageclass_model` object from a `s4t_cjs_rstan` or
`s4t_cjs_ml` object that was fit using `fixed_age = TRUE`, where the
age-class sub-model was fit separately and the age-class probabilities
were incorporated into the mark-recapture model.

## Usage

``` r
extract_ageclass_fit(object)
```

## Arguments

- object:

  a `s4t_cjs_rstan` or `s4t_cjs_ml` object

## Value

a `s4t_ageclass_model` object

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_s4t_cjs_ml(p_formula = ~ t,
                     theta_formula = ~ a1 * a2 * s * j,
                     ageclass_formula = ~ FL,
                     fixed_age = TRUE,
                     s4t_ch = sim.dat$s4t_ch)
extract_ageclass_fit(m1)
} # }
```
