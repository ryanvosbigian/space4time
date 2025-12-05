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
#> a_alpha_1             a_alpha_1  -2.7915    0.6746  -4.138 -4.1137 -1.4694
#> a_alpha_2             a_alpha_2   4.1830    0.5366   7.796  3.1313  5.2346
#> a_beta_FL             a_beta_FL   1.2206    0.1230   9.922  0.9795  1.4617
#> a_beta_obs_time a_beta_obs_time  -0.4791    0.3669  -1.306 -1.1981  0.2399
#>      AIC
#> 1 112.92
```
