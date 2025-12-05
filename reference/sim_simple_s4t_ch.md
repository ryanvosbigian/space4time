# Simulate simple capture history

Simulates a small capture history with 3 sites, where there can be
holdovers between the first 1st and 2nd sites.

## Usage

``` r
sim_simple_s4t_ch(N = 500, max_obs_year = 2, prop_missing_age = 0.5)
```

## Arguments

- N:

  An `integer` of the number of individuals

- max_obs_year:

  an `integer` of the number of total time intervals

- prop_missing_age:

  a `numeric` value between `0` and `1` of the proportion of individuals
  that do not have known ages. Default is 0.50.

## Value

a `list` containing capture history `data.frame`, `list` containing the
parameters.
