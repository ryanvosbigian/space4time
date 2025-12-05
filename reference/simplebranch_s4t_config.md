# Create site configuration for simple branching site arrangement

This is only for when there are two or more initial sites where
individuals can be released, and then the next site is for all sites.
After this shared site, there should be at least one additional site in
order to estimate transitions and apparent survival (because detection
probability cannot be estimated for the last site).

## Usage

``` r
simplebranch_s4t_config(
  sites_names,
  branch_sites,
  holdover_sites = NULL,
  min_a,
  max_a,
  sites_to_pool = NULL
)
```

## Arguments

- sites_names:

  a character `vector` of the site names. Order indicates the direction
  of movement.

- branch_sites:

  a character `vector` of the initial branching sites

- holdover_sites:

  a character `vector` of the sites where after passing, individuals can
  holdover (take 1 or more time intervals) before reaching the next
  site.

- min_a:

  a `vector` of the minimum ageclass individuals in a site can be. Must
  be the same length and order as `sites_names`

- max_a:

  a `vector` of the maximum ageclass individuals in a site can be. Must
  be the same length and order as `sites_names`

- sites_to_pool:

  a named list of character vectors that contain the names of sites to
  be pooled together and treated as one site. See LINK VIGNETTE and
  details.
