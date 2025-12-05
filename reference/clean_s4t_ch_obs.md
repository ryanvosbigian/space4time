# Clean the capture history of a space4time capture history object

Cleans the capture history of a space4time capture history object. See
details.

## Usage

``` r
clean_s4t_ch_obs(s4t_ch)
```

## Arguments

- s4t_ch:

  a `s4t_ch` object.

## Value

a `clean_s4t_ch` object.

## Details

The following issues are addressed by this function: repeat observations
of individuals at the same site, time difference in captures that
exceeds what is possible given age ranges, individuals observed after
being removed (zombies), reverse movements, and observed ages of known
age fish exceeding the minimum or maximum age for a site.

    For each of these issues, the problematic observations are removed but the
    the other observations of these individuals are retained.
