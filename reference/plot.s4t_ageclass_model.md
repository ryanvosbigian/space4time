# Plot confusion matrix for fit_ageclass objects

Plot confusion matrix for `fit_ageclass` objects of predicted age class
versus observed age class.

## Usage

``` r
# S3 method for class 's4t_ageclass_model'
plot(x, ...)
```

## Arguments

- x:

  a `s4t_ageclass_cjs` object

- ...:

  Other arguments. Not used (needed for generic consistency).

## Value

a `ggplot` figure showing the confusion matrix of the model fit.

## Examples

``` r
sim.dat <- sim_simple_s4t_ch(N = 2000)
m1 <- fit_ageclass(age_formula = ~ FL,
                   s4t_ch = sim.dat$s4t_ch)
plot(m1)



```
