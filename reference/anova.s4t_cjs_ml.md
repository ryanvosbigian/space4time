# Compute likelihood ratio tests of two fitted model objects

Compute analysis of variance table.

## Usage

``` r
# S3 method for class 's4t_cjs_ml'
anova(object, ...)
```

## Arguments

- object:

  object of class `s4t_cjs_ml`

- ...:

  object of class `s4t_cjs_ml` (for generic consistency)

## Value

Returns a data frame with the difference in degrees of freedom, test
statistic (chi-squared), and p-value for the likelihood ratio tests
between the full (assumed to be the model with the largest degrees of
freedom) and reduced model.

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
anova(m1, m2)
} # }

```
