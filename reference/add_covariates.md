# Add covariates to a space4time capture history object

Add covariates to a space4time capture history object. If the covariate
has missing levels from the indices (`a1,a2,j,k,s,t,r,g`), an indicator
variable is added that is named `"OFF_"` appended to the covariate name.

## Usage

``` r
add_covariates(cov_df, s4t_ch)
```

## Arguments

- cov_df:

  a `data.frame` or `list` of `data.frame`'s containing the covariates
  for theta and p `a1,a2,j,k,s,t,r,g` indices. See details.

- s4t_ch:

  a `s4t_ch` object.

## Value

a `s4t_ch` object with covariates added

## Details

To show how covariates are added or what levels are required for the
indices, use extract_covariates.
