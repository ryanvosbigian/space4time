
<!-- README.md is generated from README.Rmd. Please edit that file -->

# space4time

<!-- badges: start -->

[![R-CMD-check](https://github.com/ryanvosbigian/space4time/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ryanvosbigian/space4time/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of space4time is to provide a user-friendly interface for
fitting a age-specific time-stratified space-for-time mark-recapture
model.

## Installation

You can install the development version of space4time like so:

``` r
remotes::install_github("ryanvosbigian/space4time")
```

## Example

``` r
library(space4time)
## basic example code

set.seed(1)

# simulate data required to make capture history
sim.dat <- simulate_data(N = 800)

# to make the capture history, need data in a specific format.

# the observations of individuals need to be in a data.frame, with
# the columns "id", "site", "time", and "removed" (no other columns are used)
head(sim.dat$ch_df)
#> # A tibble: 6 × 4
#>      id site   time removed
#>   <dbl> <chr> <dbl> <lgl>  
#> 1     1 1         2 FALSE  
#> 2     1 2         3 FALSE  
#> 3     2 1         1 FALSE  
#> 4     2 2         1 FALSE  
#> 5     3 1         1 FALSE  
#> 6     3 2         2 FALSE
# id is the individual ID. Converted to a character type.
# site is the name of the site the individual was observed at
# time is the time stratum (i.e. year). Can be a Date format or integer.
# removed is a logical (TRUE or FALSE) that indicates where the individual
#         was removed at that site (i.e. capture but not released)

# the data on the age of the individuals 
head(sim.dat$aux_age_df)
#>      id obs_time   FL ageclass
#> [1,]  1        2 -2.8       NA
#> [2,]  2        1 -1.5        2
#> [3,]  3        1  0.9       NA
#> [4,]  4        1  4.9        3
#> [5,]  5        2 -3.6        1
#> [6,]  6        1  4.8        3
# the required columns are "id", "obs_time", and "ageclass"
# id is the individual id. Note that all individuals must be included
# obs_time is the time when the age of the individuals was observed. Or,
#           if the time when the individual was first observed.
# ageclass is the age of the individual. If age is unknown, NA is used.
# Other columns can be used in the ageclass submodel or for 
#           individual based covariates. NAs are not allowed.


# the names of the sites
print(sim.dat$sites_names)
#> [1] "1" "2" "3"

# the minimum age individuals can be at each site
print(sim.dat$min_a)
#> [1] 1 1 1

# the maximum age individuals can be at each site
print(sim.dat$max_a)
#> [1] 3 3 3

# sites_config describes the configuration of the sites. Where the row indicates
# the site an individual can move from, and column indicates the site
# an individual moves to. If individuals can move from one site to another,
# than that position in the matrix is 1, otherwise it is 0. 
# Here: individuals can move from site 1 to site 2, and from site 2 to site 3.
print(sim.dat$sites_config)
#>   1 2 3
#> 1 0 1 0
#> 2 0 0 1
#> 3 0 0 0
# The rows and columns must be in the same order of sites as sites_names. 
# Here, the names are used in the rows and column names.
# The only limitation for site configuration is that from any site, individuals
# can only move to one site. 

# sites_config describes the configuration of which sites individuals can
# holdover (wait at least one time period) before moving. It uses the same
# format sites_config, where rows indicate the site an individual is moving
# from and column is the site it is moving to. 
# Here: individuals can holdover between sites 1 and 2, but not between
# sites 2 and 3. 
print(sim.dat$holdover_config)
#>   1 2 3
#> 1 0 1 0
#> 2 0 0 0
#> 3 0 0 0


# using this information, the capture history is created:
ch <- s4t_cjs_ch(
      ch_df = sim.dat$ch_df,
      aux_age_df = sim.dat$aux_age_df,
      min_a = sim.dat$min_a,
      max_a = sim.dat$max_a,
      sites_names = sim.dat$sites_names,
      sites_config = sim.dat$sites_config,
      holdover_config = sim.dat$holdover_config
    )
#> Potential error log:
#> Number of individuals encountered more than once at a site: 0
#> Number of individuals with a gap in observation times that exceed difference in minimum and maximum ages: 0
#> Number of site/time combinations with less than 10 observations: 1
#> Number of site/time combinations with no observations: 0
#> Number of 'zombies' (individuals observed after being removed): 0
#> Number of individuals with reverse movements: 0
#> Number of known age individuals with observed ages greater than max_a: 0


# Fitting the model.
# p_formula = ~ t : detection probability varies by time (note that
#                   detection probabilities are not used for transitions from
#                   the first site or to the last site. 
#                   Because there are only 3 sites, the only site where
#                   detection probs are estimated are for site 2. 
#                   If there were additional sites, then we would use ~t*k)

# theta_formula = ~ a1 * a2 * s * j : fitting the fully saturated transition
#                                     probabilities. This allows for separate
#                                     probabilities for each age, duration of
#                                     holdover (i.e. age that individuals
#                                     transition to the next site), time, 
#                                     and site.

# fixed_age = TRUE : The age-class submodel is fit separately, and the outputs
#                     are used in the mark-recapture model. The model runs
#                     faster when this is done.

m1 <- fit_s4t_cjs_rstan(
      p_formula = ~ t,
      theta_formula = ~ a1 * a2 * s * j,
      ageclass_formula = ~ FL,
      fixed_age = TRUE,
      s4t_ch = ch,
      chains = 2,
      warmup = 200,
      iter = 600
    )
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft6d' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.001472 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 14.72 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 600 [  0%]  (Warmup)
#> Chain 1: Iteration:  60 / 600 [ 10%]  (Warmup)
#> Chain 1: Iteration: 120 / 600 [ 20%]  (Warmup)
#> Chain 1: Iteration: 180 / 600 [ 30%]  (Warmup)
#> Chain 1: Iteration: 201 / 600 [ 33%]  (Sampling)
#> Chain 1: Iteration: 260 / 600 [ 43%]  (Sampling)
#> Chain 1: Iteration: 320 / 600 [ 53%]  (Sampling)
#> Chain 1: Iteration: 380 / 600 [ 63%]  (Sampling)
#> Chain 1: Iteration: 440 / 600 [ 73%]  (Sampling)
#> Chain 1: Iteration: 500 / 600 [ 83%]  (Sampling)
#> Chain 1: Iteration: 560 / 600 [ 93%]  (Sampling)
#> Chain 1: Iteration: 600 / 600 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 33.982 seconds (Warm-up)
#> Chain 1:                45.839 seconds (Sampling)
#> Chain 1:                79.821 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft6d' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000987 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 9.87 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 600 [  0%]  (Warmup)
#> Chain 2: Iteration:  60 / 600 [ 10%]  (Warmup)
#> Chain 2: Iteration: 120 / 600 [ 20%]  (Warmup)
#> Chain 2: Iteration: 180 / 600 [ 30%]  (Warmup)
#> Chain 2: Iteration: 201 / 600 [ 33%]  (Sampling)
#> Chain 2: Iteration: 260 / 600 [ 43%]  (Sampling)
#> Chain 2: Iteration: 320 / 600 [ 53%]  (Sampling)
#> Chain 2: Iteration: 380 / 600 [ 63%]  (Sampling)
#> Chain 2: Iteration: 440 / 600 [ 73%]  (Sampling)
#> Chain 2: Iteration: 500 / 600 [ 83%]  (Sampling)
#> Chain 2: Iteration: 560 / 600 [ 93%]  (Sampling)
#> Chain 2: Iteration: 600 / 600 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 25.2 seconds (Warm-up)
#> Chain 2:                41.502 seconds (Sampling)
#> Chain 2:                66.702 seconds (Total)
#> Chain 2:


# plot the transition probabilities:
plotTheta(m1)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

# plot the apparent survival probabilities:
plotSurvival(m1)
```

<img src="man/figures/README-example-2.png" width="100%" />
