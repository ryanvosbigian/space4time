# Summarize a fitted model object

Summarize a fitted `s4t_ageclass_cjs` object.

## Usage

``` r
# S3 method for class 's4t_ageclass_model'
summary(object, ...)
```

## Arguments

- object:

  a `s4t_ageclass_cjs` object

- ...:

  Other arguments. Not used (needed for generic consistency).

## Value

A list with several fitted model quantities to create summaries when
printing, including estimated parameters AIC.

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_ageclass(ageclass_formula = ~ FL,
                   s4t_ch = sim.dat$s4t_ch)
summary(m1)
} # }

```
