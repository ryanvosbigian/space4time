# Confidence intervals for fitted model parameters for s4t_cjs_ml

Computes Wald confidence intervals for one or more parameters in a
fitted model object

## Usage

``` r
# S3 method for class 's4t_cjs_ml'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  A fitted model object of class `s4t_cjs_ml`

- parm:

  A specification of what parameters are to be given confidence
  intervals. A character vector of names.

- level:

  The confidence level. Default is `0.95`

- ...:

  additional args. Not used (needed for generic consistency)

## Value

Gaussian-based confidence intervals for the fixed effect coefficients
based on the confidence level specified by `level`. Confidence intervals
are on the logit scale.
