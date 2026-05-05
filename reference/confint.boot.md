# Bootstrap confidence intervals for Two-Level Mixed Models

This is a wrapper for getting CIs for multiple parameters after running
[`bootstrap_mer`](https://marklhc.github.io/bootmlm/reference/bootstrap_mer.md),
to be consistent with a similar method for the
[`bootMer`](https://rdrr.io/pkg/lme4/man/bootMer.html) class.

## Usage

``` r
# S3 method for class 'boot'
confint(
  object,
  parm,
  level = 0.95,
  type = c("norm", "basic", "perc", "bca"),
  L = NULL,
  ...
)
```

## Arguments

- object:

  an object returned by
  [`bootstrap_mer`](https://marklhc.github.io/bootmlm/reference/bootstrap_mer.md).

- parm:

  a specification of which parameters are to be given confidence
  intervals as a vector of numbers. If missing, all parameters are
  considered.

- level:

  the confidence level required.

- type:

  character indicating the type of intervals required, as described in
  [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html). Currently
  `"stud"` is not supported.

- L:

  empirical influence values required for `type = "bca"` as described in
  [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html).

- ...:

  additional argument(s) passed to
  [`boot.ci`](https://rdrr.io/pkg/boot/man/boot.ci.html).

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter.

## Examples

``` r
library(lme4)
fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
mySumm <- function(x) {
  c(getME(x, "beta"), sigma(x))
}

# residual bootstrap
boo_resid <- bootstrap_mer(fm01ML, mySumm, type = "residual", nsim = 100)
#> 5 occurrence(s) of: message(s) : boundary (singular) fit: see help('isSingular')
confint(boo_resid, type = "bca", L = empinf_merm(fm01ML, mySumm))
#>           2.5 %     97.5 %
#> [1,] 1497.77784 1576.23319
#> [2,]   35.37269   59.27245
```
