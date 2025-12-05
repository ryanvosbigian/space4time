# Calculate abundance estimates

Calculate abundance estimates of cohorts

## Usage

``` r
abundance_estimates(obj, abund, type = c("BroodYear", "ReleaseYear", "None"))
```

## Arguments

- obj:

  a capture history data frame (see documentation for s4t_ch)

- abund:

  a vector of sites that identify kelts (i.e. adult fish ladders)

- type:

  the type of summarization for the abundance estimates. BroodYear
  summarizes by broodyear (s - a1), ReleaseYear summarizes by time (s),
  and None does not summarize.

## Value

a capture history data frame without observations following observations
at the specified sites.
