# Extract pointwise log-likelihood from a fitted s4t_cjs_rstan model

Function for extracting and de-marginalizing the pointwise
log-likelihood matrix from a s4t_cjs_rstan object. See
loo::extract_log_lik.

## Usage

``` r
extract_log_lik_s4t(mod)
```

## Arguments

- mod:

  a `s4t_cjs_rstan` object

## Value

an `S` by `N` matrix of post-warmup extracted draws, where `S` is the
size of the posterior sample, and `N` is the number of unique
observations.
