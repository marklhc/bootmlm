# Deviance Function for Multilevel Models

Creates a deviance function for a fitted model object, using \\\theta\\
and \\sigma\\ as the parameters.

## Usage

``` r
devfun_mer(x)

devfun_mer2(x)
```

## Arguments

- x:

  A fitted merMod object from
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

## Value

A deviance function with one argument `th_sig` that takes input of a
numeric vector corresponding to the some estimated values of \\\theta\\
and \\\sigma\\. For `devfun_mer2`, a function with one argument `theta`
that profiles out \\\sigma\\ and only takes input of elements for
\\\theta\\.

## Details

The built-in function(s) in lme4 for generating deviance function are
not exported and rely on C++ code. This function is mainly used to
obtain Hessian and asymptotic covariance matrix of the random effects.

## References

Bates, D., Maechler, M., Bolker, B. M., & Walker, S. C. Fitting linear
mixed-effects models using lme4. Retrieved from
<https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf>

Implementation in lme4pureR:
<https://github.com/lme4/lme4pureR/blob/master/R/pls.R>

## Examples

``` r
library(lme4)
fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
dd <- devfun_mer(fm01ML)
# Asymptotic variance-covariance matrix of (theta, sigma):
2 * solve(numDeriv::hessian(dd, c(fm01ML@theta, sigma(fm01ML))))
#> Warning: the default value of argument 'sqrt' of method 'determinant(<CHMfactor>, <logical>)' may change from TRUE to FALSE as soon as the next release of Matrix; set 'sqrt' when programming
#>           [,1]     [,2]
#> [1,]  0.108021 -1.05037
#> [2,] -1.050370 51.06771
```
