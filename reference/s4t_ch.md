# Create space-for-time mark-recapture capture history object

Create space-for-time mark-recapture capture history object.

## Usage

``` r
s4t_ch(ch_df, aux_age_df, s4t_config, cov_p = NULL, cov_theta = NULL)
```

## Arguments

- ch_df:

  `data.frame` object containing a row for each capture event. See
  details.

- aux_age_df:

  `data.frame` containing auxiliary data for each individual. See
  details.

- s4t_config:

  a `s4t_config` object created using
  [`s4t_config()`](https://ryanvosbigian.github.io/space4time/reference/s4t_config.md),
  `linear_s4t_config`, or `simplebranch_s4t_config`.

- cov_p:

  a `data.frame` or `list` of `data.frame`'s containing the covariates
  for p `a1,a2,j,k,s,t,r,g` indices. See details.

- cov_theta:

  a `data.frame` or `list` of `data.frame`'s containing the covariates
  for theta `a1,a2,j,k,s,t,r,g` indices. See details.

## Details

The capture history data (`ch_df`) must be a `data.frame` (or coercible
to a `data.frame`) with exactly four columns named `id`, `site`, `time`,
and `removed`. Each row is a release and recapture (or observation)
event. The `id` column is the unique identifier for each individual. The
`site` column is the site name, which must correspond to the names in
the `s4t_config` object. The `time` column must either be an integer for
the time period or a Date. If it is a date, it will be converted to
years. The last column is a `logical`(i.e. `TRUE` or `FALSE`) that
indicates whether individuals were removed (i.e. retained) at the event.

The auxiliary and age data (`aux_age_df`) must be a `data.frame` (or
coercible to a `data.frame`) that must contain at least three columns
named `id`, `obs_time`, and `ageclass`. Additional columns can be
included that contain data on individuals. The `id` column is the unique
identifier for individuals, `obs_time` is the integer time period (or
Date) when the individual was first observed or when the `ageclass` of
the individual was observed. `ageclass` is the integer age of the
individual. If the `ageclass` was not observed, then `ageclass = NA`,
but `obs_time` must be filled in. `obs_time` should correspond to the
time period of the auxiliary data.

The `cov_p` and `cov_theta` arguments can be used to add in covariates.
They can be data frames or a list of data frames. The data frames are
joined (using
[`dplyr::left_join`](https://dplyr.tidyverse.org/reference/mutate-joins.html))
to the \$theta\$ or \$p\$ parameters using indices for site, time,
initial release group, age, and group. Better practice is to use the
[`add_covariates()`](https://ryanvosbigian.github.io/space4time/reference/add_covariates.md)
function to add any covariates rather than adding it manually, so that
any missing levels can be addressed. To see the indices of the
parameters, use
[`extract_covariates()`](https://ryanvosbigian.github.io/space4time/reference/extract_covariates.md).

Note that individual covariates can be included in the `s4t_cjs_ml` and
`s4t_cjs_rstan` models. These covariates are included in the
`aux_age_df` data.

## Examples

``` r
ch_df <- data.frame(id = c(1,1,1,
                           2,2,
                           3,3,
                           4,
                           5,
                           6),
                     site = c("A","B","C",
                              "A","B",
                              "A","C",
                              "A",
                              "A",
                              "A"),
                     time = c(1,3,3,
                              2,3,
                              1,3,
                              2,
                              1,
                              1),
                     removed = c(FALSE,FALSE,FALSE,
                                 FALSE,FALSE,
                                 FALSE,FALSE,
                                 FALSE,
                                 FALSE,
                                 FALSE)
                      )

aux_age_df <- data.frame(id = 1:6,
                          obs_site = rep("A",6),
                          ageclass = c(1,2,1,1,2,1),
                          obs_time = c(1,2,1,2,1,1),
                          Covariate1 = c(3,1,2,1,2,1))

site_arr <- linear_s4t_config(sites_names = c("A","B","C"),
                              holdover_sites = c("A"),
                              min_a = c(1,1,1),
                              max_a = c(3,3,3))

ch <- s4t_ch(ch_df = ch_df,
             aux_age_df = aux_age_df,
             s4t_config = site_arr)
#> 
#> Error log:
#> 
#> Repeat encounters at same site N = 0
#> Individuals observed after being removed ('zombies') N = 0
#> Gap in observation times that exceed max difference in ages N = 0
#> Holdovers observed between sites with only direct transitions N = 0
#> Reverse movements N = 0
#> Known age individuals with ages outside of site-specific age-range N = 0
#> Individuals with missing initial release site N = 0
#> 
#> Potential errors:
#> Site/time combinations with 0 observations N = 0
#> Site/time combinations with less than 10 observations N = 4



```
