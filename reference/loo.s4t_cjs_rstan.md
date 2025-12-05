# Efficient approximate leave-one-out cross-validation (LOO)

Computes the PSIS-LOO CV efficient approximate leave-one-out (LOO)
cross-validation for Bayesian models using Pareto smoothed importance
sampling (PSIS). See documentation for
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)

## Usage

``` r
# S3 method for class 's4t_cjs_rstan'
loo(x, ...)
```

## Arguments

- x:

  a `s4t_cjs_rstan` object

- ...:

  passed to loo::loo

## Value

a named list with class `c("psis_loo", "loo")`. See
[`?loo::loo`](https://mc-stan.org/loo/reference/loo.html)
