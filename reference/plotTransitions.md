# Plot transition probability estimates from fitted `s4t_cjs_ml` or `s4t_cjs_rstan`objects

Plot transition probability estimates (which include movement and
survival probabilities) from fitted `s4t_cjs_ml` or
`s4t_cjs_rstan`objects. If the package `geomtextpath` is installed, the
figure is more visually appealing.

## Usage

``` r
plotTransitions(x, textsize = 3, ...)
```

## Arguments

- x:

  a `s4t_cjs` or `s4t_cjs_rstan` object

- textsize:

  an integer for the font size to pass to
  [`geomtextpath::geom_textsegment()`](https://allancameron.github.io/geomtextpath/reference/geom_textsegment.html)

- ...:

  passed to
  [`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html)
  for selecting particular survivals. Filters apparent survivals from
  `s4t_cjs_ml$apparent_surv` or `s4t_cjs_rstan$apparent_surv`

## Value

a `ggplot2` figure showing estimated transition probabilities.

## Examples

``` r
if (FALSE) { # \dontrun{
 sim.dat <- sim_simple_s4t_ch(N = 2000)
 m1 <- fit_s4t_cjs_rstan(p_formula = ~ t,
                         theta_formula = ~ a1 * a2 * s * j,
                         ageclass_formula = ageclass ~ FL,
                         fixed_age = TRUE,
                         s4t_ch = sim.dat$s4t_ch)
plotTransitions(m1)
} # }
```
