# Markov chain traceplots

Draw the traceplot corresponding to one or more Markov chains, providing
a visual way to inspect sampling behavior and assess mixing across
chains and convergence.

## Usage

``` r
traceplot.s4t_cjs_rstan(object, pars = NULL, ...)
```

## Arguments

- object:

  a `s4t_cjs_rstan` object

- pars:

  a `character` object that parses parameters using regular expressions.
  Defaults to all parameters or the first 10 parameters ( if there are
  more than 10)

- ...:

  Optional arguments to pass to `stan::traceplot()`

## Value

A `ggplot` object that can be further customized using the `ggplot2`
package.
