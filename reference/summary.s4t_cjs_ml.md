# Summarize a fitted model object

Summarize a fitted `s4t_cjs_ml` object.

## Usage

``` r
# S3 method for class 's4t_cjs_ml'
summary(object, ...)
```

## Arguments

- object:

  a `s4t_cjs_ml` object

- ...:

  Other arguments. Not used (needed for generic consistency).

## Value

A list with several fitted model quantities to create summaries when
printing, including estimated parameters and AIC.

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_s4t_cjs_ml(p_formula = ~ t,
                      theta_formula = ~ a1 * a2 * s * j,
                      ageclass_formula = ~ FL,
                      fixed_age = TRUE,
                      s4t_ch = sim.dat$s4t_ch)
summary(m1)
} # }
```
