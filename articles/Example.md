# Example using PTAGIS queries

``` r
library(space4time)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

This is an example workflow for implementing the model. This
implementation uses data queried from Columbia Basin Research Data
Access in Real Time (DART). The Columbia Basin Research has a query that
was developed for use with the Basin TribPit software. However, the data
can also be used by *space4time*.

We start with the observations of individuals from the East Fork
Potlatch River Rotary Screw Trap, which has the site name: EFPTRP
(previously POTREF). We’ll use observations from 2015 to 2019 for this
example. These data are from Idaho Department of Fish and Game’s
juvenile trapping database.

Here are some of the rows of a file that contain data on juvenile
steelhead encountered at this rotary screw trap. Some of the columns
(`id`, `obs_time`, and `ageclass`) are required to have those names. Any
additional variables can be included. However, complete data for all
variables (except `ageclass`) is required for any of the following
analyses.

``` r
EFPTRP_age_data[10:15,]%>%
  knitr::kable(row.names = FALSE)
```

| id             | obs_time | SurveyDate |  FL | Bin | Season | ageclass |
|:---------------|---------:|:-----------|----:|----:|:-------|---------:|
| 3DD.007753945C |     2015 | 5/3/2015   |  80 |  80 | Spring |       NA |
| 3DD.007753732D |     2015 | 5/4/2015   |  80 |  80 | Spring |        1 |
| 3DD.0077535753 |     2015 | 5/5/2015   |  80 |  80 | Spring |       NA |
| 3DD.0077532003 |     2015 | 5/8/2015   |  80 |  80 | Spring |       NA |
| 3DD.0077537A2A |     2015 | 5/9/2015   |  80 |  80 | Spring |       NA |
| 3DD.007752C3C2 |     2015 | 5/11/2015  |  80 |  80 | Spring |       NA |

We identify and drop duplicate observations as well as drop individuals
with missing auxilary data (`FL`).

``` r

# identify individuals with more than one observation at the RST in the data
EFPTRP_age_data %>%
  dplyr::group_by(id) %>% # groups data by ID
  dplyr::arrange(id,SurveyDate) %>% # sorts by ID and date
  dplyr::mutate(Obs = dplyr::n()) %>% # creates a column that is the number of times
                                      # each individual is in the data
  filter(Obs > 1) # filter for multiple observations
#> # A tibble: 0 × 8
#> # Groups:   id [0]
#> # ℹ 8 variables: id <chr>, obs_time <int>, SurveyDate <chr>, FL <int>,
#> #   Bin <int>, Season <chr>, ageclass <int>, Obs <int>

# identify missing data
EFPTRP_age_data %>%
  dplyr::filter(is.na(FL)) # filter for FL as NA
#> [1] id         obs_time   SurveyDate FL         Bin        Season     ageclass  
#> <0 rows> (or 0-length row.names)

# drop duplicate observations (keep the first observation) and missing data
EFPTRP_age_data <- EFPTRP_age_data %>%
  dplyr::group_by(id) %>%
  dplyr::arrange(id,SurveyDate) %>%
  dplyr::summarize(         # summarizes by id, so it returns 1 row per id
            id = dplyr::first(id), # retains the first value
            SurveyDate = dplyr::first(SurveyDate),
            obs_time = dplyr::first(obs_time), 
            FL = dplyr::first(FL),
            Bin = dplyr::first(Bin),
            Season = dplyr::first(Season),
            ageclass = dplyr::first(ageclass)
            # obs_site = first(obs_site), # save the first observation of any additional variables
            ) %>%
  dplyr::filter(!is.na(FL)) 
```

Instead of using fork lengths, for this example we will use binned fork
lengths. This allows for greater efficency in model fitting because it
allows for observations to be marginalized (grouped together). The bins
are also z-scored (centered and divided by the standard deviation) to
help with model fitting if Bin is treated as a continuous variable. The
min and max bin are truncated so that there are not any bins with only a
handful of observations.

``` r
# create histogram to visualize the distribution of fork lengths
# hist(EFPTRP_age_data$FL,breaks = 50)


EFPTRP_age_data <- EFPTRP_age_data %>%
  dplyr::mutate(Bin = dplyr::case_when(FL < 80 ~ 80, # if FL < 80, return 80
                                       FL >180 ~ 180, # if FL greater than 180, return 180
                                       # otherwise, round FL to the nearest 10 mm.
                                       FL >= 80 & FL <= 180 ~ round(FL,-1)), 
                Bin_sc = (Bin - mean(Bin)) / stats::sd(Bin)) # z-scale Bin
```

To obtain capture history data, we used query tools developed by DART
for Basin TribPit using functions available in the `space4time` package.
We used tools developed for Basin TribPit and PitPro to query PTAGIS and
identify transported fish. More information regarding Basin TribPit and
the associated queries can be found at
“<https://www.cbr.washington.edu/analysis/apps/BasinTribPit>”. To
conduct queries, we used the query “Upload TagID List created by the
user” (“<https://www.cbr.washington.edu/dart/query/pit_tagids>”). We
uploaded the tag list of observations from the age data and conducted
the query to create a “TribPit Observation File”
(“pitbasin_branching_upload_1759171262_907.csv”). To identify
transported fish, we use their “site_config.txt” file, which located on
the page for PitPro
(“<https://www.cbr.washington.edu/analysis/apps/pitpro>”) and is
regularly updated
(“<https://www.cbr.washington.edu/paramest/docs/pitpro/updates/sites_config.txt>”).

An alternative to using the queries in DART is to use a Complete Tag
History query in PTAGIS. However, care must be taken to identify
indiviudals that were transported by barges.

The function to process DART TribPit observation files is
[`read_DART_file()`](https://ryanvosbigian.github.io/space4time/reference/read_DART_file.md).
The first argument is the filepath to the “TribPit Observation File”,
the second is the age data, that need to the specific columns identified
above, and third (optional) is the path to the “site_config.txt” file.
The default is to use the link.

``` r
proc_DART_data <- read_DART_file(filepath = "pitbasin_branching_upload_1759171262_907.csv",
               aux_age_df = EFPTRP_age_data,
               DART_config = "sites_config.txt")
#> Obtained the PTAGIS site configuration formatted by DART from:
#> sites_config.txt
#> Warning in identify_barged_fish(capture_data, parsed_df): Number of
#> observations without info. Assuming not transported. N = 1
```

Note see
[`?read_DART_file`](https://ryanvosbigian.github.io/space4time/reference/read_DART_file.md)
if there were more than one DART file to include.

The output from `read_DART_file` (`proc_DART_data`) is a list with two
objects, a data.frame of capture occasions (`proc_DART_data$ch_df`) and
a data frame of the age auxiliary data (`proc_DART_data$aux_age_df`)

``` r
ch_df <- proc_DART_data$ch_df
head(ch_df)
#>               id   site time removed
#> 1 384.3B23958521    BCC 2015   FALSE
#> 2 384.3B23958521 POTREF 2015   FALSE
#> 3 384.3B23958521    TWX 2015   FALSE
#> 4 384.3B2397A68E    GOJ 2015   FALSE
#> 5 384.3B2397A68E    GRJ 2015   FALSE
#> 6 384.3B2397A68E POTREF 2015   FALSE

aux_age_df <- proc_DART_data$aux_age_df

head(aux_age_df)
#> # A tibble: 6 × 8
#>   id             SurveyDate obs_time    FL   Bin Season ageclass Bin_sc
#>   <chr>          <chr>         <int> <int> <dbl> <chr>     <int>  <dbl>
#> 1 384.3B23958521 3/15/2015      2015   205   180 Spring       NA   2.47
#> 2 384.3B2397A68E 3/19/2015      2015   160   160 Spring       NA   1.74
#> 3 384.3B239C9845 3/20/2015      2015   177   180 Spring       NA   2.47
#> 4 384.3B239D8C19 3/19/2015      2015   160   160 Spring       NA   1.74
#> 5 384.3B239EB35C 3/19/2015      2015   137   140 Spring        2   1.01
#> 6 384.3B239ED4C7 3/19/2015      2015   138   140 Spring       NA   1.01
```

Summarize some of the observations

``` r
table(aux_age_df$ageclass)
#> 
#>   1   2   3 
#> 388 192  21
```

``` r
table(ch_df$site)
#> 
#>    B2J    BBA    BCC    BO1    BO2    BO3    BO4    DRM EFPTRP    EPR    GOA 
#>      3      1     38      2      2      3      5      1    685      1      5 
#>    GOJ    GRA    GRJ    GRS    HLM    ICH    JDJ    JO1    JO2    KHS    LAP 
#>    126      5    191      8    259     21     39      2      1      5      1 
#>    LMA    LMJ    MC1    MCJ    PCM POTREF    TD1    TD2    TWX 
#>      5     84      5     19      2   2451      4      2     14
```

In this case, some of the observations are likely of kelts, which are
post-spawning adult steelhead attempting to return to the ocean. To
exclude these, we can use the
[`remove_kelt_obs()`](https://ryanvosbigian.github.io/space4time/reference/remove_kelt_obs.md)
function.

The `kelt_obssite` argument is the site after which all observations
should be dropped, meaning that if an individual is observed at this
site, then all further observations of the individual in the same or
next time periods should be dropped. The Lower Granite Dam fish ladder
site (“GRA”) is used to identify adults.

``` r
ch_df2 <- remove_kelt_obs(ch_df = ch_df,kelt_obssite = "GRA")
#> Dropping observations after observations at the following sites: GRA
```

Next, we can drop observations at non-target sites.

``` r
# only keep observations in the following sites
ch_df3 <- ch_df2 %>% 
  dplyr::filter(site %in% c("EFPTRP","POTREF","GOJ","GRJ","LMJ","MCJ","TWX","JDJ","BCC","B2J")) 
```

Note that here, the name of the RST in the EF Potlatch River was POTREF
at some points in time.

Next, we can check what age range to include for each site:

``` r
# returns a data frame of the ages of known age fish at different observations
ch_df %>%
  dplyr::left_join(aux_age_df,by = "id") %>%  # merge age data with capture history
  dplyr::filter(!is.na(ageclass)) %>% # only check observations of known age fish
  dplyr::mutate(Age_at_obs = ageclass + (time - obs_time)) %>% # calculate known age
  dplyr::group_by(site,Age_at_obs) %>% # group by site and age at site
  dplyr::summarise(N = dplyr::n()) # summarize number of observations
#> # A tibble: 32 × 3
#> # Groups:   site [19]
#>    site   Age_at_obs     N
#>    <chr>       <dbl> <int>
#>  1 BCC             2     2
#>  2 BO2             3     1
#>  3 BO4             3     1
#>  4 EFPTRP          1    76
#>  5 EFPTRP          2    59
#>  6 EFPTRP          3     7
#>  7 GOA             3     1
#>  8 GOJ             1     3
#>  9 GOJ             2    25
#> 10 GOJ             3     3
#> # ℹ 22 more rows
```

The range of ages at each site is 1 through 3. All sites with
substantial observations have the same range.

This ends the initial data cleaning. Now, we can create site
configuration object. Because there are multiple names for the initial
release site, these are treated as if there are multiple release sites
so we use the `simplebranch_s4t_config` function. The sites are
reprinted below.

``` r
table(ch_df3$site)
#> 
#>    B2J    BCC EFPTRP    GOJ    GRJ    JDJ    LMJ    MCJ POTREF    TWX 
#>      3     38    685    126    191     39     84     19   2451     14
```

See documentation for more details
([`?simplebranch_s4t_config`](https://ryanvosbigian.github.io/space4time/reference/simplebranch_s4t_config.md)).
The East Fork Potlatch RST has been at multiple locations, although here
we treat them as the same site because they are not at substantially
different locations. The sites_to_pool argument is used to merge these
sites together. The holdover sites are the sites after which individuals
can holdover before transitioning to the next site. The `branch_sites`
argument specifies the initial branching sites, which are the two names
for the East Fork Potlatch RST.

Here, we assume that no fish holdover after passing Lower Granite Dam
(“GRJ”). We pool all recapture sites below Little Goose Dam (“GOJ”) with
Little Goose. We have to specify min age and max age for each site even
if they are pooled.

``` r
ef_pot_site_config <- simplebranch_s4t_config(sites_names = c("EFPTRP",
                                                        "POTREF",
                                                        "GRJ","GOJ","LMJ","MCJ","JDJ",
                                                        "BCC",
                                                        "TWX"),
                                              branch_sites = c("EFPTRP",
                                                          "POTREF"),
                                       holdover_sites = c("EFPTRP",
                                                          "POTREF"),
                                       sites_to_pool = list("EFPTRP" = c("EFPTRP","POTREF"),
                                                            "GOJ" = c("LMJ","MCJ","JDJ",
                                                                      "BCC","TWX")),
                                       min_a = c("EFPTRP" = 1,
                                                 "POTREF" = 1,
                                                 "GRJ" = 1,
                                                 "GOJ" = 1,
                                                 "LMJ" = 1,
                                                 "MCJ" = 1,
                                                 "JDJ" = 1,
                                                 "BCC" = 1,
                                                 "TWX" = 1),
                                       max_a = c("EFPTRP" = 3,
                                                 "POTREF" = 3,
                                                 "GRJ" = 3,
                                                 "GOJ" = 3,
                                                 "LMJ" = 3,
                                                 "MCJ" = 3,
                                                 "JDJ" = 3,
                                                 "BCC" = 3,
                                                 "TWX" = 3)

)
```

Next, we create initial capture history object:

``` r
efp_s4_ch <- s4t_ch(ch_df = ch_df3,aux_age_df = aux_age_df,s4t_config = ef_pot_site_config)
#> Removing sites from capture history: B2J
#> Note, there are IDs in aux_age_df that are not in ch_df, n = 1
#> 
#> Error log:
#> 
#> Repeat encounters at same site N = 0
#> Individuals observed after being removed ('zombies') N = 3
#> Gap in observation times that exceed max difference in ages N = 3
#> Holdovers observed between sites with only direct transitions N = 0
#> Reverse movements N = 1
#> Known age individuals with ages outside of site-specific age-range N = 0
#> Individuals with missing initial release site N = 0
#> 
#> Potential errors:
#> Site/time combinations with 0 observations N = 0
#> Site/time combinations with less than 10 observations N = 2
#> Maximum release occasions in l_matrix and m_matrix do not match.
```

We need to address the errors in the error log. We can use the
[`clean_s4t_ch_obs()`](https://ryanvosbigian.github.io/space4time/reference/clean_s4t_ch_obs.md)
function, which returns a cleaned capture history data frame as well as
information on what observations were dropped.

``` r
clean_ch_df <- clean_s4t_ch_obs(efp_s4_ch)
#> N = 4 observations and N = 0 individuals were dropped.
```

We can inspect what observations were dropped and determine whether
these makes biological sense.

``` r

# show dropped observations
clean_ch_df$dropped_ch_df

# show full observation history of individuals with dropped observations
clean_ch_df$intermediate_ch_df %>%
  dplyr::group_by(id) %>% # group by ID
  dplyr::filter(any(drop_obs == TRUE)) %>% # keep only individuals with dropped observations
  head() # print out only the first 6 observations


# can inspect individuals

# aux_age_df %>%
#   dplyr::filter(id == "3DD.00778C962F")

# clean_ch_df$intermediate_ch_df %>%
#   dplyr::filter(id == "3DD.00778C962F")
```

Next, we re-make the capture history object using the cleaned `ch_df`
data frame.

``` r
efp_s4_ch2 <- s4t_ch(ch_df = clean_ch_df$cleaned_ch_df,aux_age_df = aux_age_df,s4t_config = ef_pot_site_config)
#> Note, there are IDs in aux_age_df that are not in ch_df, n = 1
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
#> Site/time combinations with less than 10 observations N = 1
```

There are no errors, so we can use this capture history object in our
analyses. If there were more errors (which does occasionally happen), we
could clean the `efp_s4_ch2` object and see if the errors can be
corrected by the cleaning function.

Next, we can fit the model for age class. First, we will compute some
summaries of the ageing data.

``` r
# annual summaries
aux_age_df %>%
  dplyr::filter(!is.na(ageclass)) %>% # only include observations with ages
  dplyr::group_by(obs_time,ageclass) %>% # group by age and obs_time (year)
  dplyr::summarize(mean_FL = mean(FL), # mean FL
            sd_FL = stats::sd(FL),   # sd of FL
            N = dplyr::n())          # number of observations
#> `summarise()` has grouped output by 'obs_time'. You can override using the
#> `.groups` argument.
#> # A tibble: 14 × 5
#> # Groups:   obs_time [5]
#>    obs_time ageclass mean_FL sd_FL     N
#>       <int>    <int>   <dbl> <dbl> <int>
#>  1     2015        1    98.0 12.4     91
#>  2     2015        2   146.  16.4     53
#>  3     2015        3   176.  15.9      9
#>  4     2016        1   103.  12.2     57
#>  5     2016        2   144.  16.2     41
#>  6     2016        3   176.   7.22     5
#>  7     2017        1   106.  15.2    164
#>  8     2017        2   150.  19.2     39
#>  9     2018        1    99.2 12.4     51
#> 10     2018        2   146.  23.2     31
#> 11     2018        3   148   NA        1
#> 12     2019        1    99.8 12.6     25
#> 13     2019        2   148.  16.0     28
#> 14     2019        3   166.  12.0      6


# create boxplot of fork length by ages
aux_age_df %>% 
  dplyr::filter(!is.na(ageclass)) %>% # only include observations with ages
  ggplot2::ggplot(ggplot2::aes(factor(ageclass), FL)) +
  geom_boxplot() +
  facet_wrap(~obs_time) +
  theme_bw() +
  labs(x = "Age",y = "Fork length (mm)")
```

![](Example_files/figure-html/unnamed-chunk-19-1.png)

Based on the above figure, we have the following models for ageclass,
which use ordinal regression. We could use FL or scaled FL, but we will
instead use the binned fork length and the scaled binned fork length.
This is so that the model can be more efficient as a result of
marginalization.

``` r

age_mod1 <- fit_ageclass(age_formula = ~ I(factor(obs_time)) + Bin_sc,
                         s4t_ch = efp_s4_ch2)

age_mod2 <- fit_ageclass(age_formula = ~ I(factor(obs_time)) * Bin_sc,
                         s4t_ch = efp_s4_ch2)

age_mod3 <- fit_ageclass(age_formula = ~ I(factor(obs_time)) + I(factor(Bin)),
                         s4t_ch = efp_s4_ch2)
```

We can inspect some simple goodness of fit:

``` r
plot(age_mod1)
```

![](Example_files/figure-html/unnamed-chunk-21-1.png)

If further goodness-of-fit metrics are desired, there is a `predict`
function for the ageclass models that can be used to return predicted
values.

We can use AIC to select the best fitting model:

``` r
AIC(age_mod1,age_mod2,age_mod3)
#>          df      AIC
#> age_mod1  7 177.8796
#> age_mod2 11 181.8134
#> age_mod3 16 190.2178
```

The first model is selected based on AIC.

Next, we fit space-for-time mark recapture models. The ageclass model is
the top model from above (note: use one-sided formulas only). The
formulas for detection probability (`p`) and conditional transition
probability (`theta`) are determined based on the site structure. There
is only one site, so the full model for theta is
`theta ~ a1 * a2 * s * j`, which allows for different transitions based
on time at “release” (`s`), “release” site (`j`), age at “release”
(`a1`), and age at “recapture” (`a2`). The quotes around “release” and
“recapture” indicate that the captures can be passive and that in this
case does not imply actual capture (i.e. if they passed a site but
weren’t detected). We recommend fitting the full model for theta to
obtain transition estimates for every combination of age, time, and
site. Many reduced formulas would not make sense.

The full model for p is `p ~ t * a1 * a2`, which allows for different
detection probability for each time of “recapture” (`t`), age at
“release” (`a1`), and age at “recapture” (`a2`). Site is technically
included, but because detection probabilities are not estimated for the
first or final site, only the second site (“GRJ”) has detection
probabilities estimated. The detection probabilities for the final site
are not separable from transition rates (note: this means that
transition rates estimated between “GRJ” and “GOJ” are not actually the
transition rates).

The full model for p may be over-parameterized, so we also fit reduced
formulas where detection probability only depends on time and where
detection probability depends on time and age at “recapture”.

``` r

s4t_m1 <- fit_s4t_cjs_rstan(p_formula = ~ t,
                            theta_formula = ~ a1 * a2 * s * j,
                            ageclass_formula = ~ I(factor(obs_time)) + Bin_sc,
                            s4t_ch = efp_s4_ch2)
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000822 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 8.22 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 70.818 seconds (Warm-up)
#> Chain 1:                77.32 seconds (Sampling)
#> Chain 1:                148.138 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000653 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 6.53 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 60.846 seconds (Warm-up)
#> Chain 2:                39.346 seconds (Sampling)
#> Chain 2:                100.192 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000723 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 7.23 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 78.47 seconds (Warm-up)
#> Chain 3:                56.989 seconds (Sampling)
#> Chain 3:                135.459 seconds (Total)
#> Chain 3:


s4t_m2 <- fit_s4t_cjs_rstan(p_formula = ~ t * a2,
                            theta_formula = ~ a1 * a2 * s * j,
                            ageclass_formula = ~ I(factor(obs_time)) + Bin_sc,
                            s4t_ch = efp_s4_ch2)
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000826 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 8.26 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 81.699 seconds (Warm-up)
#> Chain 1:                57.77 seconds (Sampling)
#> Chain 1:                139.469 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000688 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 6.88 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 89.663 seconds (Warm-up)
#> Chain 2:                58.536 seconds (Sampling)
#> Chain 2:                148.199 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.000683 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 6.83 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 80.824 seconds (Warm-up)
#> Chain 3:                58.237 seconds (Sampling)
#> Chain 3:                139.061 seconds (Total)
#> Chain 3:

# full model
s4t_m3 <- fit_s4t_cjs_rstan(p_formula = ~ t * a1 * a2,
                            theta_formula = ~ a1 * a2 * s * j,
                            ageclass_formula = ~ I(factor(obs_time)) + Bin_sc,
                            s4t_ch = efp_s4_ch2)
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.0015 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 15 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 82.091 seconds (Warm-up)
#> Chain 1:                57.999 seconds (Sampling)
#> Chain 1:                140.09 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0.000677 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 6.77 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 2: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 2: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 2: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 85.531 seconds (Warm-up)
#> Chain 2:                60.896 seconds (Sampling)
#> Chain 2:                146.427 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 's4t_cjs_fixedage_draft7' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0.001356 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 13.56 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 3: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 3: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 3: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 3: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 96.527 seconds (Warm-up)
#> Chain 3:                66.899 seconds (Sampling)
#> Chain 3:                163.426 seconds (Total)
#> Chain 3:

# optionally save model objects
# save(s4t_m1,s4t_m2,s4t_m3, file = "EF_Potlatch_s4t.RData")
```

We can check model fits:

``` r
s4t_m1
#> Inference for Stan model: s4t_cjs_fixedage_draft7.
#> 3 chains, each with iter=1000; warmup=500; thin=1; 
#> post-warmup draws per chain=500, total post-warmup draws=1500.
#> 
#>                       mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff
#> theta_(Intercept)    -4.56    0.02 0.49 -5.59 -4.87 -4.56 -4.22 -3.64   603
#> theta_a12             0.77    0.02 0.78 -0.76  0.24  0.76  1.29  2.32  1291
#> theta_a13             4.30    0.02 0.82  2.61  3.79  4.32  4.85  5.85  1405
#> theta_a22             1.90    0.02 0.51  0.92  1.55  1.90  2.22  2.94   705
#> theta_a23            -1.07    0.02 0.63 -2.38 -1.48 -1.06 -0.66  0.15   844
#> theta_s2016           0.51    0.02 0.70 -1.00  0.07  0.53  0.99  1.79   854
#> theta_s2017           0.59    0.02 0.55 -0.51  0.23  0.58  0.95  1.70   673
#> theta_s2018          -0.33    0.02 0.77 -1.89 -0.83 -0.30  0.18  1.09  1126
#> theta_s2019          -0.16    0.02 0.85 -1.88 -0.72 -0.13  0.40  1.49  1190
#> theta_s2020           1.03    0.03 1.16 -1.18  0.27  0.99  1.82  3.34  1836
#> theta_jGRJ            0.40    0.02 0.81 -1.21 -0.15  0.42  0.96  1.92  1267
#> theta_a12:a22         2.02    0.02 0.77  0.48  1.51  2.02  2.52  3.60  1042
#> theta_a12:s2016       0.09    0.03 1.16 -2.20 -0.70  0.10  0.89  2.38  1548
#> theta_a13:s2016       0.61    0.03 1.02 -1.38 -0.10  0.61  1.28  2.61  1568
#> theta_a12:s2017       1.27    0.02 0.98 -0.57  0.58  1.26  1.94  3.24  1662
#> theta_a13:s2017       0.12    0.03 1.45 -2.85 -0.80  0.10  1.03  3.24  2947
#> theta_a12:s2018      -0.03    0.03 1.09 -2.07 -0.78 -0.08  0.69  2.18  1675
#> theta_a13:s2018       1.00    0.03 1.25 -1.59  0.23  1.05  1.88  3.27  1416
#> theta_a12:s2019      -0.53    0.03 0.93 -2.43 -1.11 -0.53  0.07  1.34   940
#> theta_a13:s2019       0.33    0.03 1.16 -2.06 -0.44  0.33  1.13  2.62  1631
#> theta_a12:s2020       0.85    0.03 1.26 -1.44  0.02  0.78  1.67  3.39  1786
#> theta_a22:s2016       0.71    0.03 0.76 -0.71  0.20  0.69  1.20  2.27   875
#> theta_a23:s2016       0.55    0.03 0.88 -1.18 -0.02  0.54  1.12  2.31   895
#> theta_a22:s2017      -0.16    0.02 0.62 -1.42 -0.56 -0.15  0.26  1.01   727
#> theta_a23:s2017       0.32    0.03 0.76 -1.14 -0.20  0.28  0.86  1.92   914
#> theta_a22:s2018       0.96    0.02 0.80 -0.56  0.42  0.95  1.49  2.60  1142
#> theta_a23:s2018       0.56    0.03 1.02 -1.43 -0.12  0.59  1.28  2.45  1416
#> theta_a22:s2019       0.99    0.03 0.96 -0.78  0.32  0.93  1.64  2.89  1205
#> theta_a23:s2019      -0.01    0.03 1.12 -2.20 -0.77  0.01  0.75  2.18  1565
#> theta_a12:jGRJ       -1.98    0.03 0.96 -3.88 -2.56 -1.99 -1.32 -0.06   835
#> theta_a13:jGRJ        0.27    0.03 1.10 -1.84 -0.49  0.26  1.01  2.44  1752
#> theta_s2016:jGRJ     -0.83    0.02 1.05 -2.96 -1.58 -0.82 -0.09  1.19  1800
#> theta_s2017:jGRJ      1.18    0.03 0.92 -0.64  0.57  1.17  1.76  3.03  1112
#> theta_s2018:jGRJ     -0.22    0.03 1.13 -2.46 -0.98 -0.22  0.52  2.07  1621
#> theta_s2019:jGRJ      0.77    0.03 1.07 -1.36  0.05  0.78  1.45  2.86  1218
#> theta_a12:a22:s2016  -0.01    0.03 1.10 -2.00 -0.76 -0.04  0.68  2.39  1461
#> theta_a12:a22:s2017   0.40    0.02 0.92 -1.40 -0.26  0.43  1.06  2.11  1478
#> theta_a12:a22:s2018   0.56    0.03 1.09 -1.49 -0.18  0.52  1.32  2.74  1766
#> theta_a12:s2016:jGRJ -0.63    0.03 1.18 -2.91 -1.45 -0.64  0.18  1.69  1644
#> theta_a13:s2016:jGRJ  0.11    0.03 1.32 -2.37 -0.78  0.14  1.02  2.63  2468
#> theta_a12:s2017:jGRJ -1.79    0.03 1.06 -3.85 -2.47 -1.83 -1.06  0.37   976
#> theta_a13:s2017:jGRJ  0.05    0.03 1.41 -2.72 -0.85  0.04  1.00  2.71  2385
#> theta_a12:s2018:jGRJ  0.09    0.03 1.20 -2.32 -0.74  0.10  0.88  2.47  1727
#> theta_a13:s2018:jGRJ -0.20    0.03 1.28 -2.65 -1.07 -0.22  0.66  2.37  2096
#> theta_a12:s2019:jGRJ  0.09    0.03 1.20 -2.37 -0.67  0.10  0.90  2.48  1308
#> theta_a13:s2019:jGRJ  0.94    0.03 1.33 -1.57  0.06  0.89  1.84  3.55  1922
#> p_(Intercept)        -1.89    0.02 0.37 -2.61 -2.15 -1.89 -1.65 -1.16   434
#> p_t2016               0.03    0.02 0.44 -0.85 -0.28  0.04  0.34  0.86   628
#> p_t2017               0.95    0.02 0.40  0.18  0.66  0.96  1.22  1.70   451
#> p_t2018               0.83    0.02 0.43 -0.01  0.54  0.84  1.13  1.69   636
#> p_t2019               1.15    0.02 0.45  0.29  0.86  1.15  1.45  2.03   560
#> p_t2020              -3.85    0.06 2.26 -9.19 -5.14 -3.59 -2.24 -0.38  1267
#> overall_surv[1]       0.08    0.00 0.03  0.04  0.06  0.08  0.10  0.14  1688
#> overall_surv[2]       0.53    0.01 0.13  0.31  0.43  0.52  0.62  0.81   572
#> overall_surv[3]       0.23    0.00 0.12  0.05  0.15  0.21  0.30  0.50  1749
#> overall_surv[4]       0.22    0.00 0.04  0.16  0.20  0.22  0.24  0.30  1552
#> overall_surv[5]       0.77    0.00 0.15  0.45  0.67  0.79  0.89  0.99  1332
#> overall_surv[6]       0.56    0.00 0.20  0.19  0.41  0.55  0.72  0.93  1920
#> overall_surv[7]       0.13    0.00 0.02  0.09  0.11  0.12  0.14  0.18  1689
#> overall_surv[8]       0.90    0.00 0.06  0.76  0.86  0.90  0.94  0.99  2144
#> overall_surv[9]       0.45    0.01 0.29  0.03  0.19  0.43  0.70  0.96  2066
#> overall_surv[10]      0.13    0.00 0.04  0.08  0.11  0.13  0.15  0.22  1628
#> overall_surv[11]      0.75    0.00 0.14  0.49  0.64  0.75  0.86  0.98  1684
#> overall_surv[12]      0.50    0.01 0.27  0.02  0.30  0.54  0.72  0.93  1360
#> overall_surv[13]      0.17    0.00 0.11  0.05  0.10  0.14  0.20  0.49  1176
#> overall_surv[14]      0.59    0.00 0.14  0.37  0.49  0.57  0.67  0.93  1223
#> overall_surv[15]      0.27    0.00 0.17  0.03  0.15  0.25  0.38  0.66  1692
#> overall_surv[16]      0.02    0.00 0.02  0.00  0.01  0.02  0.03  0.07  1282
#> overall_surv[17]      0.20    0.00 0.06  0.11  0.15  0.19  0.23  0.33   775
#> overall_surv[18]      0.36    0.00 0.17  0.09  0.23  0.34  0.47  0.73  1701
#> overall_surv[19]      0.02    0.00 0.03  0.00  0.00  0.01  0.03  0.09  1920
#> overall_surv[20]      0.17    0.00 0.05  0.10  0.14  0.17  0.20  0.31  1481
#> overall_surv[21]      0.54    0.01 0.22  0.17  0.37  0.52  0.71  0.96  1779
#> overall_surv[22]      0.10    0.00 0.06  0.02  0.05  0.09  0.13  0.26  2020
#> overall_surv[23]      0.51    0.00 0.06  0.41  0.47  0.51  0.55  0.64  2116
#> overall_surv[24]      0.71    0.01 0.30  0.05  0.51  0.83  0.96  1.00  1354
#> overall_surv[25]      0.02    0.00 0.03  0.00  0.00  0.01  0.02  0.11  1331
#> overall_surv[26]      0.40    0.00 0.10  0.24  0.33  0.39  0.46  0.63  1882
#> overall_surv[27]      0.53    0.01 0.26  0.05  0.33  0.55  0.75  0.95  1407
#> overall_surv[28]      0.06    0.00 0.08  0.00  0.01  0.03  0.07  0.28  1480
#> overall_surv[29]      0.43    0.00 0.11  0.24  0.35  0.43  0.50  0.67  1530
#> overall_surv[30]      0.72    0.00 0.18  0.32  0.59  0.74  0.87  0.98  1709
#> overall_surv[31]      0.57    0.01 0.24  0.14  0.38  0.57  0.78  0.97  1647
#> overall_surv[32]      0.56    0.01 0.25  0.11  0.35  0.58  0.78  0.96  1674
#> cohort_surv[1]        0.01    0.00 0.01  0.00  0.01  0.01  0.01  0.03   676
#> cohort_surv[2]        0.07    0.00 0.02  0.03  0.05  0.06  0.08  0.12  1923
#> cohort_surv[3]        0.00    0.00 0.00  0.00  0.00  0.00  0.00  0.01  1434
#> cohort_surv[4]        0.53    0.01 0.13  0.31  0.43  0.52  0.62  0.80   571
#> cohort_surv[5]        0.00    0.00 0.00  0.00  0.00  0.00  0.01  0.02  1534
#> cohort_surv[6]        0.23    0.00 0.12  0.05  0.15  0.21  0.30  0.50  1749
#> cohort_surv[7]        0.02    0.00 0.02  0.00  0.01  0.02  0.03  0.06  1642
#> cohort_surv[8]        0.19    0.00 0.03  0.14  0.17  0.19  0.21  0.25  1598
#> cohort_surv[9]        0.01    0.00 0.01  0.00  0.01  0.01  0.01  0.03  1809
#> cohort_surv[10]       0.76    0.00 0.15  0.44  0.66  0.78  0.89  0.98  1295
#> cohort_surv[11]       0.01    0.00 0.01  0.00  0.00  0.00  0.01  0.03  1980
#> cohort_surv[12]       0.56    0.00 0.20  0.19  0.41  0.55  0.72  0.93  1920
#> cohort_surv[13]       0.02    0.00 0.01  0.01  0.01  0.02  0.02  0.04  2180
#> cohort_surv[14]       0.10    0.00 0.02  0.06  0.08  0.09  0.11  0.15  1472
#> cohort_surv[15]       0.01    0.00 0.00  0.00  0.01  0.01  0.01  0.02  1699
#> cohort_surv[16]       0.89    0.00 0.06  0.75  0.85  0.89  0.93  0.98  2187
#> cohort_surv[17]       0.01    0.00 0.01  0.00  0.00  0.01  0.01  0.03  1772
#> cohort_surv[18]       0.45    0.01 0.29  0.03  0.19  0.43  0.70  0.96  2066
#> cohort_surv[19]       0.01    0.00 0.01  0.00  0.00  0.01  0.01  0.03  1437
#> cohort_surv[20]       0.12    0.00 0.03  0.07  0.09  0.11  0.14  0.19  1544
#> cohort_surv[21]       0.01    0.00 0.01  0.00  0.00  0.00  0.01  0.02  1933
#> cohort_surv[22]       0.74    0.00 0.14  0.48  0.64  0.74  0.86  0.98  1682
#> cohort_surv[23]       0.00    0.00 0.01  0.00  0.00  0.00  0.00  0.02  1670
#> cohort_surv[24]       0.50    0.01 0.27  0.02  0.30  0.54  0.72  0.93  1360
#> cohort_surv[25]       0.01    0.00 0.01  0.00  0.00  0.01  0.02  0.04  1555
#> cohort_surv[26]       0.16    0.00 0.11  0.05  0.09  0.12  0.19  0.47  1135
#> cohort_surv[27]       0.59    0.00 0.14  0.37  0.49  0.56  0.67  0.92  1225
#> cohort_surv[28]       0.00    0.00 0.00  0.00  0.00  0.00  0.00  0.02  1510
#> cohort_surv[29]       0.27    0.00 0.17  0.03  0.15  0.25  0.38  0.66  1692
#> cohort_surv[30]       0.02    0.00 0.02  0.00  0.01  0.02  0.03  0.07  1282
#> cohort_surv[31]       0.20    0.00 0.06  0.11  0.15  0.19  0.23  0.33   775
#> cohort_surv[32]       0.36    0.00 0.17  0.09  0.23  0.34  0.47  0.73  1701
#> cohort_surv[33]       0.02    0.00 0.03  0.00  0.00  0.01  0.03  0.09  1920
#> cohort_surv[34]       0.17    0.00 0.05  0.10  0.14  0.17  0.20  0.31  1481
#> cohort_surv[35]       0.54    0.01 0.22  0.17  0.37  0.52  0.71  0.96  1779
#> cohort_surv[36]       0.10    0.00 0.06  0.02  0.05  0.09  0.13  0.26  2020
#> cohort_surv[37]       0.51    0.00 0.06  0.41  0.47  0.51  0.55  0.64  2116
#> cohort_surv[38]       0.71    0.01 0.30  0.05  0.51  0.83  0.96  1.00  1354
#> cohort_surv[39]       0.02    0.00 0.03  0.00  0.00  0.01  0.02  0.11  1331
#> cohort_surv[40]       0.40    0.00 0.10  0.24  0.33  0.39  0.46  0.63  1882
#> cohort_surv[41]       0.53    0.01 0.26  0.05  0.33  0.55  0.75  0.95  1407
#> cohort_surv[42]       0.06    0.00 0.08  0.00  0.01  0.03  0.07  0.28  1480
#> cohort_surv[43]       0.43    0.00 0.11  0.24  0.35  0.43  0.50  0.67  1530
#> cohort_surv[44]       0.72    0.00 0.18  0.32  0.59  0.74  0.87  0.98  1709
#> cohort_surv[45]       0.57    0.01 0.24  0.14  0.38  0.57  0.78  0.97  1647
#> cohort_surv[46]       0.56    0.01 0.25  0.11  0.35  0.58  0.78  0.96  1674
#>                      Rhat
#> theta_(Intercept)    1.00
#> theta_a12            1.00
#> theta_a13            1.00
#> theta_a22            1.00
#> theta_a23            1.00
#> theta_s2016          1.00
#> theta_s2017          1.00
#> theta_s2018          1.00
#> theta_s2019          1.00
#> theta_s2020          1.00
#> theta_jGRJ           1.00
#> theta_a12:a22        1.00
#> theta_a12:s2016      1.00
#> theta_a13:s2016      1.00
#> theta_a12:s2017      1.00
#> theta_a13:s2017      1.00
#> theta_a12:s2018      1.00
#> theta_a13:s2018      1.00
#> theta_a12:s2019      1.00
#> theta_a13:s2019      1.00
#> theta_a12:s2020      1.00
#> theta_a22:s2016      1.00
#> theta_a23:s2016      1.00
#> theta_a22:s2017      1.00
#> theta_a23:s2017      1.00
#> theta_a22:s2018      1.00
#> theta_a23:s2018      1.00
#> theta_a22:s2019      1.00
#> theta_a23:s2019      1.00
#> theta_a12:jGRJ       1.00
#> theta_a13:jGRJ       1.00
#> theta_s2016:jGRJ     1.00
#> theta_s2017:jGRJ     1.00
#> theta_s2018:jGRJ     1.00
#> theta_s2019:jGRJ     1.00
#> theta_a12:a22:s2016  1.00
#> theta_a12:a22:s2017  1.00
#> theta_a12:a22:s2018  1.00
#> theta_a12:s2016:jGRJ 1.00
#> theta_a13:s2016:jGRJ 1.00
#> theta_a12:s2017:jGRJ 1.01
#> theta_a13:s2017:jGRJ 1.00
#> theta_a12:s2018:jGRJ 1.00
#> theta_a13:s2018:jGRJ 1.00
#> theta_a12:s2019:jGRJ 1.00
#> theta_a13:s2019:jGRJ 1.00
#> p_(Intercept)        1.01
#> p_t2016              1.00
#> p_t2017              1.01
#> p_t2018              1.01
#> p_t2019              1.00
#> p_t2020              1.00
#> overall_surv[1]      1.00
#> overall_surv[2]      1.01
#> overall_surv[3]      1.00
#> overall_surv[4]      1.00
#> overall_surv[5]      1.00
#> overall_surv[6]      1.00
#> overall_surv[7]      1.00
#> overall_surv[8]      1.00
#> overall_surv[9]      1.00
#> overall_surv[10]     1.00
#> overall_surv[11]     1.00
#> overall_surv[12]     1.00
#> overall_surv[13]     1.00
#> overall_surv[14]     1.00
#> overall_surv[15]     1.00
#> overall_surv[16]     1.00
#> overall_surv[17]     1.00
#> overall_surv[18]     1.00
#> overall_surv[19]     1.00
#> overall_surv[20]     1.00
#> overall_surv[21]     1.00
#> overall_surv[22]     1.00
#> overall_surv[23]     1.00
#> overall_surv[24]     1.00
#> overall_surv[25]     1.00
#> overall_surv[26]     1.00
#> overall_surv[27]     1.00
#> overall_surv[28]     1.00
#> overall_surv[29]     1.00
#> overall_surv[30]     1.00
#> overall_surv[31]     1.00
#> overall_surv[32]     1.00
#> cohort_surv[1]       1.00
#> cohort_surv[2]       1.00
#> cohort_surv[3]       1.00
#> cohort_surv[4]       1.01
#> cohort_surv[5]       1.00
#> cohort_surv[6]       1.00
#> cohort_surv[7]       1.00
#> cohort_surv[8]       1.00
#> cohort_surv[9]       1.00
#> cohort_surv[10]      1.00
#> cohort_surv[11]      1.00
#> cohort_surv[12]      1.00
#> cohort_surv[13]      1.00
#> cohort_surv[14]      1.00
#> cohort_surv[15]      1.00
#> cohort_surv[16]      1.00
#> cohort_surv[17]      1.00
#> cohort_surv[18]      1.00
#> cohort_surv[19]      1.00
#> cohort_surv[20]      1.00
#> cohort_surv[21]      1.00
#> cohort_surv[22]      1.00
#> cohort_surv[23]      1.00
#> cohort_surv[24]      1.00
#> cohort_surv[25]      1.00
#> cohort_surv[26]      1.00
#> cohort_surv[27]      1.00
#> cohort_surv[28]      1.00
#> cohort_surv[29]      1.00
#> cohort_surv[30]      1.00
#> cohort_surv[31]      1.00
#> cohort_surv[32]      1.00
#> cohort_surv[33]      1.00
#> cohort_surv[34]      1.00
#> cohort_surv[35]      1.00
#> cohort_surv[36]      1.00
#> cohort_surv[37]      1.00
#> cohort_surv[38]      1.00
#> cohort_surv[39]      1.00
#> cohort_surv[40]      1.00
#> cohort_surv[41]      1.00
#> cohort_surv[42]      1.00
#> cohort_surv[43]      1.00
#> cohort_surv[44]      1.00
#> cohort_surv[45]      1.00
#> cohort_surv[46]      1.00
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Dec  4 16:47:29 2025.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Note that the R-hats for all parameters should all be below 1.05. It is
also important to inspect traceplots to check that the chains are
properly mixed:

``` r
traceplot.s4t_cjs_rstan(s4t_m1)
```

![](Example_files/figure-html/unnamed-chunk-25-1.png)

We can compare the models using loo-psis:

``` r
library(loo)
#> This is loo version 2.8.0
#> - Online documentation and vignettes at mc-stan.org/loo
#> - As of v2.0.0 loo defaults to 1 core but we recommend using as many as possible. Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session.
#> - Windows 10 users: loo may be very slow if 'mc.cores' is set in your .Rprofile file (see https://github.com/stan-dev/loo/issues/94).


loo_m1 <- loo(s4t_m1)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
loo_m2 <- loo(s4t_m2)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
loo_m3 <- loo(s4t_m3)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

loo::loo_compare(loo_m1,loo_m2,loo_m3)
#>        elpd_diff se_diff
#> model1   0.0       0.0  
#> model2  -0.6       5.7  
#> model3 -10.5       7.2
```

The best model is the first model, although the first and second are
similar. The third is also a reasonable model, but it has slightly worse
fit.

We can extract the cohort transition rates, which are the probability
that an individual transitions at a particular time.

``` r
s4t_m1$cohort_transitions %>% head()
#>                a1 a2    s    t      j   k      r g        mean      se_mean
#> cohort_surv[1]  1  1 2015 2015 EFPTRP GRJ EFPTRP 1 0.011596141 2.255094e-04
#> cohort_surv[2]  1  2 2015 2016 EFPTRP GRJ EFPTRP 1 0.067762938 5.202243e-04
#> cohort_surv[3]  1  3 2015 2017 EFPTRP GRJ EFPTRP 1 0.003859254 6.247052e-05
#> cohort_surv[4]  2  2 2015 2015 EFPTRP GRJ EFPTRP 1 0.527944021 5.450687e-03
#> cohort_surv[5]  2  3 2015 2016 EFPTRP GRJ EFPTRP 1 0.004703763 1.013366e-04
#> cohort_surv[6]  3  3 2015 2015 EFPTRP GRJ EFPTRP 1 0.230996056 2.796183e-03
#>                         sd         2.5%         25%         50%         75%
#> cohort_surv[1] 0.005861567 0.0037382181 0.007585499 0.010358318 0.014463397
#> cohort_surv[2] 0.022811372 0.0326458232 0.051413544 0.064563167 0.081115477
#> cohort_surv[3] 0.002365250 0.0010690766 0.002282107 0.003337353 0.004890911
#> cohort_surv[4] 0.130286186 0.3061972524 0.429073455 0.516816675 0.615960422
#> cohort_surv[5] 0.003969004 0.0006248236 0.002107785 0.003656135 0.006084381
#> cohort_surv[6] 0.116954616 0.0535492085 0.145051610 0.214602857 0.302718925
#>                      97.5%     n_eff      Rhat
#> cohort_surv[1] 0.025672911  675.6133 1.0043057
#> cohort_surv[2] 0.120210019 1922.7446 1.0012341
#> cohort_surv[3] 0.009525344 1433.5208 0.9998680
#> cohort_surv[4] 0.803037193  571.3395 1.0066268
#> cohort_surv[5] 0.015675809 1534.0167 0.9993725
#> cohort_surv[6] 0.496569670 1749.4581 0.9986282
```

We can also extract the apparent survival, which are the probability
that an individual transitions at any time (sum of the cohort_transition
rates for a particular age, time, and site):

``` r
s4t_m1$apparent_surv %>% head()
#>                 a1    s      j   k      r g       mean      se_mean         sd
#> overall_surv[1]  1 2015 EFPTRP GRJ EFPTRP 1 0.08321833 0.0006211894 0.02552496
#> overall_surv[2]  2 2015 EFPTRP GRJ EFPTRP 1 0.53264778 0.0054519093 0.13038041
#> overall_surv[3]  3 2015 EFPTRP GRJ EFPTRP 1 0.23099606 0.0027961830 0.11695462
#> overall_surv[4]  1 2016 EFPTRP GRJ EFPTRP 1 0.22043826 0.0009113620 0.03590648
#> overall_surv[5]  2 2016 EFPTRP GRJ EFPTRP 1 0.77100705 0.0040430397 0.14755549
#> overall_surv[6]  3 2016 EFPTRP GRJ EFPTRP 1 0.56280593 0.0046261495 0.20269724
#>                       2.5%        25%        50%        75%     97.5%     n_eff
#> overall_surv[1] 0.04380003 0.06441708 0.08005725 0.09814767 0.1417122 1688.4264
#> overall_surv[2] 0.31356379 0.43454555 0.52254028 0.61892467 0.8079113  571.9097
#> overall_surv[3] 0.05354921 0.14505161 0.21460286 0.30271893 0.4965697 1749.4581
#> overall_surv[4] 0.15739545 0.19654483 0.21712918 0.24104679 0.2996030 1552.2579
#> overall_surv[5] 0.44954688 0.66714630 0.78910629 0.89423648 0.9856612 1331.9709
#> overall_surv[6] 0.18673929 0.40968979 0.55348538 0.71659939 0.9312145 1919.8016
#>                      Rhat
#> overall_surv[1] 1.0023496
#> overall_surv[2] 1.0066533
#> overall_surv[3] 0.9986282
#> overall_surv[4] 1.0017769
#> overall_surv[5] 1.0004279
#> overall_surv[6] 0.9982456
```

We can visualize these transition rates:

``` r
# currently the textsize must be set to include filters, which are j == "EFPTRP"
# in this example
plotTransitions(s4t_m1,textsize = 3,j == "EFPTRP")
#> Warning: There were 2 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `site_diff = as.integer(as.character(k)) -
#>   as.integer(as.character(j))`.
#> Caused by warning:
#> ! NAs introduced by coercion
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.
```

![](Example_files/figure-html/unnamed-chunk-29-1.png)

We can also visaulize the apparent survival rates:

``` r
plotSurvival(s4t_m1,textsize = 3,j == "EFPTRP")
#> Warning: There were 2 warnings in `dplyr::mutate()`.
#> The first warning was:
#> ℹ In argument: `site_diff = as.integer(as.character(k)) -
#>   as.integer(as.character(j))`.
#> Caused by warning:
#> ! NAs introduced by coercion
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 1 remaining warning.
```

![](Example_files/figure-html/unnamed-chunk-30-1.png)

The helper function
[`abundance_estimates()`](https://ryanvosbigian.github.io/space4time/reference/abundance_estimates.md)
can be used to calculate abundance of individuals that transition from
one site to another site. An example of where this would be useful is if
there are age-specific abundance estimates at East Fork Potlatch RST for
each year. The abundances (and standard errors of the abundance
estimates if available) are specified in the argument `abund` and the
summarization is specified by `type`. The format of the `abund` data
frame is shown below, where it needs columns for site (`j`), time (`s`),
age at release (`a1`), and abundance. This is merged with the
`cohort_transitions` data frame and the abundances are estimated.
Additionally, standard errors are approximated assuming that the
covariances of the cohort transitions are multivariate normally
distributed and that abundance estimates are independent.

``` r
# format of abund argument using fake (simulated) data.
head(fake_abundance_data)
#>        j    s a1 abundance abundance_se
#> 1 EFPTRP 2015  1       471            0
#> 2 EFPTRP 2015  2       355            0
#> 3 EFPTRP 2015  3       384            0
#> 4 EFPTRP 2016  1       157            0
#> 5 EFPTRP 2016  2       282            0
#> 6 EFPTRP 2016  3       773            0

# compute cohort specific abundances
cohort_abundance_at_LGR = abundance_estimates(s4t_m1,abund = fake_abundance_data,type = "None")
#> Setting Nan values in vcov matrix to 0, results are approximate.

# summarize by broodyear, site, and group.
broodyear_abundance_at_LGR = abundance_estimates(s4t_m1,abund = fake_abundance_data,type = "BroodYear")
#> Setting Nan values in vcov matrix to 0, results are approximate.

head(broodyear_abundance_at_LGR)
#> # A tibble: 6 × 6
#> # Groups:   broodyear, r, g [3]
#>   broodyear r          g j      abundance_broodyear abundance_broodyear_se
#>       <dbl> <chr>  <dbl> <chr>                <dbl>                  <dbl>
#> 1      2012 EFPTRP     1 EFPTRP                88.7                   44.9
#> 2      2012 EFPTRP     1 GRJ                   NA                      0  
#> 3      2013 EFPTRP     1 EFPTRP               624.                   165. 
#> 4      2013 EFPTRP     1 GRJ                   NA                      0  
#> 5      2014 EFPTRP     1 EFPTRP               627.                   239. 
#> 6      2014 EFPTRP     1 GRJ                   NA                      0
```
