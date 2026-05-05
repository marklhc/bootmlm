# Profile Likelihood Confidence Interval for Intraclass Correlation

Compute confidence intervals for the intraclass correlation of a model
fit of class
[`merMod-class`](https://rdrr.io/pkg/lme4/man/merMod-class.html).

## Usage

``` r
prof_ci_icc(x, level = 0.95, eps_max = 1e+06)
```

## Arguments

- x:

  A fitted merMod object from
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

- level:

  Confidence level between 0 and 1. Default is .95.

- eps_max:

  The maximum value that the upper limit can be. Default is 1e6.

## Value

A named numeric vector with two values showing the lower and upper limit
of the intraclass correlation.

## Details

This function uses the [`uniroot`](https://rdrr.io/r/stats/uniroot.html)
function and determine the lower and upper limit by evaluating the
profile deviance (for ML) or the -2 profile REML criterion. It works by
obtaining the interval for \\\theta = \tau / \sigma\\ and transforming
the two limits.

The resulting interval is bounded by zero for its lower limit.
