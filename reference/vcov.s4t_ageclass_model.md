# Calculate variance-covariance matrix for a fitted model object

Calculate variance-covariance matrix for a `s4t_ageclass_model` model
object.

## Usage

``` r
# S3 method for class 's4t_ageclass_model'
vcov(object, ...)
```

## Arguments

- object:

  a `s4t_ageclass_model` object fitted using
  [`fit_ageclass()`](https://ryanvosbigian.github.io/space4time/reference/fit_ageclass.md).

- ...:

  Other arguments. Not used (needed for generic consistency)

## Value

The variance-covariance matrix of coefficients from a fitted
`s4t_ageclass_model` object

## Examples

``` r
if (FALSE) { # \dontrun{
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_ageclass(ageclass_formula = ~ FL,
                   s4t_ch = sim.dat$s4t_ch)

vcov(m1)
} # }
```
