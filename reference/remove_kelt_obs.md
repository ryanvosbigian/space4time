# Remove observations of kelts

Removes any observations following an observation at specified sites.
Note that the observations at the specified sites are retained.

## Usage

``` r
remove_kelt_obs(ch_df, kelt_obssite)
```

## Arguments

- ch_df:

  a capture history data frame (see documentation for s4t_ch)

- kelt_obssite:

  a vector of sites that identify kelts (i.e. adult fish ladders)

## Value

a capture history data frame without observations following observations
at the specified sites.
