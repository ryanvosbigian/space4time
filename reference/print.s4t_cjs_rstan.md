# Print summary of `s4t_cjs_rstan`

Print summary of `s4t_cjs_rstan`

## Usage

``` r
# S3 method for class 's4t_cjs_rstan'
print(
  x,
  pars = NULL,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  digits_summary = 2,
  include = TRUE,
  ...
)
```

## Arguments

- x:

  s4t_cjs_rstan object

- pars:

  vector of parameters to include (defaults to all parameters)

- probs:

  quantiles for credible intervals

- digits_summary:

  the number of significant digits to use when printing the summary

- include:

  logical scalar indicating whether to include the parameters named by
  the pars argument.

- ...:

  passed to `summary` method for s4t_cjs_rstan object.
