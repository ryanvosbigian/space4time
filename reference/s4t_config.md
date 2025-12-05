# Create space4time site configuration object (`s4t_config`)

Create the space4time site configuration object (`s4t_config`) using
custom site configuration. For a simple configuration where there is
only one initial release site use
[`linear_s4t_config()`](https://ryanvosbigian.github.io/space4time/reference/linear_s4t_config.md),
which assumes that the sites are arranged in a simple sequence. If there
are multiple initial release sites that are immediately connected to a
simple sequence then
[`simplebranch_s4t_config()`](https://ryanvosbigian.github.io/space4time/reference/simplebranch_s4t_config.md)
can be used. Otherwise, use this function to create a custom site
configuration. See details or (LINK VIGNETTE) for more detailed
instructions.

## Usage

``` r
s4t_config(
  sites_names,
  sites_config,
  holdover_config,
  min_a,
  max_a,
  sites_to_pool = NULL
)
```

## Arguments

- sites_names:

  a character `vector` of the site names

- sites_config:

  a `matrix` that determine how sites are linked together. Must be in
  the same order as sites_names. See LINK VIGNETTE and details.

- holdover_config:

  a `matrix` that determine whether individuals can holdover in between
  sites. In the same format as sites_configs. See LINK VIGNETTE and
  details.

- min_a:

  a `vector` of the minimum ageclass individuals in a site can be. Must
  be the same length and order as `sites_names`

- max_a:

  a `vector` of the maximum ageclass individuals in a site can be. Must
  be the same length and order as `sites_names`

- sites_to_pool:

  a named list of character vectors that contain the names of sites to
  be pooled together and treated as one site. See examples, LINK
  VIGNETTE, and details.

## Value

A `s4t_config` object, which contains

## Details

Additional details... `sites_to_pool` is either a named list or a list
of character vectors. If it is unnamed, than the name, min_a, max_a
values are used from the first site listed. If the list is named, than
the name, min_a, and max_a are used from the site in the name.

## Examples

``` r
site_arrangement <- s4t_config(sites_names = c("A","B","C","D"),
                               min_a = c(1,1,1,1),
                               max_a = c(3,4,4,4),
                               sites_config = matrix(
                                       c(0,1,0,0,
                                         0,0,1,0,
                                         0,0,0,1,
                                         0,0,0,0),
                                         nrow = 4,
                                         ncol = 4,
                                         byrow = TRUE,
                                         dimnames = list(c("A","B","C","D"),
                                         c("A","B","C","D"))),
                               holdover_config = matrix(
                                        c(0,1,0,0,
                                          0,0,0,0,
                                          0,0,0,0,
                                          0,0,0,0),
                                          nrow = 4,
                                          ncol = 4,
                                          byrow = TRUE,
                                          dimnames = list(c("A","B","C","D"),
                                          c("A","B","C","D")))
)

```
