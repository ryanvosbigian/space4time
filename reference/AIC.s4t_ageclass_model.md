# Compute and return AIC of fitted model objects

Compute AIC of one or more fitted `s4t_ageclass_model` model objects

## Usage

``` r
# S3 method for class 's4t_ageclass_model'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  A `s4t_ageclass_model` model object

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
m1 <- fit_ageclass(ageclass_formula = ~ FL + obs_time,
                   s4t_ch = sim.dat$s4t_ch)
m2 <- fit_ageclass(ageclass_formula = ~ FL + I(FL^2) + obs_time,
                   s4t_ch = sim.dat$s4t_ch)
AIC(m1, m2)
} # }
```
