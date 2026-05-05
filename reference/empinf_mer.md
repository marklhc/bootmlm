# Empirical Influence Values for Two-Level Mixed Models

This function calculates the empirical influence values for a statistic
in a given fitted model object using the delete-\\m_j\\ jackknife.

## Usage

``` r
empinf_mer(x, FUN, index = 1)

empinf_merm(x, FUN)
```

## Arguments

- x:

  A fitted merMod object from
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

- FUN:

  A function taking a fitted merMod object as input and returning the
  statistic of interest.

- index:

  An integer stating the position of the statistic in the output of
  `FUN(x)`.

## Value

A numeric vector with length equals to number of clusters of `x`
containing the weighted influence value of each cluster.

## Details

`empinf_mer` computes non-parametric influence function of models fitted
using [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) by deleting one
cluster at a time. See van der Leeden, Meijer, and Busing (2008, pp.
420–422) for more information. Whereas `empinf_mer` computes influence
values for a specified position (as specified with the `index` argument)
of the output of `FUN`, `empinf_merm` computes influence values for
every element in `FUN(x)`.

## References

Van der Leeden, R., Meijer, E., & Busing, F. M. T. A. (2008). Resampling
multilevel models. In J. de Leeuw & E. Meijer (Eds.), Handbook of
multilevel Analysis (pp. 401–433). New York, NY: Springer.

## Examples

``` r
library(lme4)
fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
# Define function for intraclass correlation
icc <- function(x) 1 / (1 + 1 / getME(x, "theta")^2)
empinf_mer(fm01ML, icc)
#> [1] -2.5346963 -1.1239525 -0.1245598 -2.7668690  4.2457632  2.3043145
empinf_mer(fm01ML, fixef)
#> [1] -112.5    2.5  182.5 -147.5  362.5 -287.5
```
