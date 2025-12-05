# Compute likelihood ratio tests of two fitted model objects

Compute analysis of variance table.

## Usage

``` r
# S3 method for class 's4t_ageclass_model'
anova(object, ...)
```

## Arguments

- object:

  object of class `s4t_ageclass_model` created using
  [`fit_ageclass()`](https://ryanvosbigian.github.io/space4time/reference/fit_ageclass.md)

- ...:

  object of class `s4t_ageclass_model` (for generic consistency)

## Value

Returns a data frame with the difference in degrees of freedom, test
statistic (chi-squared), and p-value for the likelihood ratio tests
between the full (assumed to be the model with the largest degrees of
freedom) and reduced model.

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_ageclass(ageclass_formula = ~ FL + obs_time,
                   s4t_ch = sim.dat$s4t_ch)
m2 <- fit_ageclass(ageclass_formula = ~ FL + I(FL^2) + obs_time,
                   s4t_ch = sim.dat$s4t_ch)
anova(m1, m2)
} # }

```
