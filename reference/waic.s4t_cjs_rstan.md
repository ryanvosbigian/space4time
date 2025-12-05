# Widely applicable information criterion (WAIC)

The waic() methods can be used to compute WAIC from the pointwise
log-likelihood. See
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html)

## Usage

``` r
# S3 method for class 's4t_cjs_rstan'
waic(x, ...)
```

## Arguments

- x:

  a `s4t_cjs_rstan` object

- ...:

  not used (for generic consistency)

## Value

A named list (of class `c("waic", "loo")`). See
[`?loo::waic`](https://mc-stan.org/loo/reference/waic.html).
