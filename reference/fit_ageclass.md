# Fit age-class model using ordinal regression

Fit age-class ordinal regression model with logit link and flexible
threshold parameters.

## Usage

``` r
fit_ageclass(age_formula = ~1, s4t_ch)
```

## Arguments

- age_formula:

  a one- or two-sided `formula` describing the effect structure of the
  model. Note that obs_time is treated as factor

- s4t_ch:

  a `s4t_ch` object

## Value

a fitted `s4t_ageclass_model` object that contains the estimated
parameters, AIC, negative log likelihood, and call.

## Details

To enforce obs_time to be treated as a continuous variable, the formula
can be written as: `age_formula = ~ I(numeric(obs_time))`, or another
variable for time can be included in auxiliary data.

## Examples

``` r
sim.dat <- sim_simple_s4t_ch()
m <- fit_ageclass(ageclass ~ FL + obs_time,s4t_ch = sim.dat$s4t_ch)
summary(m)
#>                       parameter estimate std_error z_value   lcl95   ucl95
#> a_alpha_1             a_alpha_1   -3.604    0.6641  -5.427 -4.9056 -2.3022
#> a_alpha_2             a_alpha_2    4.198    0.5031   8.345  3.2121  5.1843
#> a_beta_FL             a_beta_FL    1.136    0.1095  10.377  0.9214  1.3505
#> a_beta_obs_time a_beta_obs_time   -1.050    0.3588  -2.926 -1.7529 -0.3465
#>      AIC
#> 1 122.81
```
