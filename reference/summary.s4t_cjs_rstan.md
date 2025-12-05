# Summary method for `s4t_cjs_rstan` objects

Summarize the distributions of the estimated parameters and derived
quantities regarding MCMC sampling. Calls `summary,stanfit-method`.

## Usage

``` r
# S3 method for class 's4t_cjs_rstan'
summary(object, pars = NULL, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
```

## Arguments

- object:

  A `s4t_cjs_rstan` object.

- pars:

  A `character` value that returns parameters that match using regular
  expressions. Defaults to all parameters.

- probs:

  a numerical vector of quantiles of interest. Default is
  `c(0.025, 0.25, 0.50, 0.75, 0.975)`

- ...:

  Additional arguments to pass to `summary,stanfit-method`.

## Value

Returns a named list with elements `summary` and `c_summary`, which
contain the summaries for all chains merged and individual chains,
respectively. See `summary,stanfit-method` documentation for more
details.

## Examples

``` r
if (FALSE) { # \dontrun{
 sim.dat <- sim_simple_s4t_ch(N = 2000)
 m1 <- fit_s4t_cjs_rstan(p_formula = ~ t,
                         theta_formula = ~ a1 * a2 * s * j,
                         ageclass_formula = ageclass ~ FL,
                         fixed_age = TRUE,
                         s4t_ch = sim.dat$s4t_ch)
summary(m1)
} # }

```
