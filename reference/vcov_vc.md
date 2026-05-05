# Asymptotic Covariance Matrix for Random Effects

Return the asymptotic covariance matrix of random effect standard
deviations (or variances) for a fitted model object, using the Hessian
evaluated at the (restricted) maximum likelihood estimates.

## Usage

``` r
vcov_vc(x, sd_cor = TRUE, print_names = TRUE)
```

## Arguments

- x:

  A fitted merMod object from
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

- sd_cor:

  Logical indicating whether to return asymptotic covariance matrix on
  SD scale (if `TRUE`) or on variance scale (if `FALSE`).

- print_names:

  Logical, whether to print the names for the covariance matrix.

## Value

A (q + 1) \* (q + 1) symmetric matrix of the covariance matrix of
(\\\tau, \sigma\\) (if `sd_cor = TRUE`) or (\\\tau^2, \sigma^2\\) (if
`sd_cor = FALSE`), where q is the the number of estimated random-effects
components (excluding \\\sigma\\). For example, for a model with random
slope, \\\tau\\ = (intercept SD, intercept-slope correlation, slope SD).

## Details

Although it's easy to obtain the Hessian for \\\theta\\, the relative
Cholesky factor, in lme4, there is no easy way to obtain the Hessian for
the variance components. This function uses
[`devfun_mer()`](https://marklhc.github.io/bootmlm/reference/devfun_mer.md)
to obtain the Hessian (\\H\\) of variance components (or standard
deviations, SD), and then obtain the asymptotic covariance matrix as
\\-2 H^{-1}\\.

## See also

[`vcov.merMod`](https://rdrr.io/pkg/lme4/man/vcov.merMod.html) for
covariance matrix of fixed effects,
[`confint.merMod`](https://rdrr.io/pkg/lme4/man/confint.merMod.html) for
confidence intervals of all parameter estimates, and
[`devfun_mer`](https://marklhc.github.io/bootmlm/reference/devfun_mer.md)
for the underlying function to produce the deviance function.

## Examples

``` r
library(lme4)
data(Orthodont, package = "nlme")
fm1 <- lmer(distance ~ age + (age | Subject), data = Orthodont)
vc <- VarCorr(fm1)
# Standard deviation only
print(vc, comp = c("Std.Dev"))
#>  Groups   Name        Std.Dev. Corr   
#>  Subject  (Intercept) 2.32736         
#>           age         0.22645  -0.609 
#>  Residual             1.31002         
# Asymptotic variance-covariance matrix of (tau, sigma):
vcov_vc(fm1, sd_cor = TRUE)
#>                             sd_(Intercept)|Subject cor_age.(Intercept)|Subject
#> sd_(Intercept)|Subject                  1.13460585                 -0.29955866
#> cor_age.(Intercept)|Subject            -0.29955866                  0.10590733
#> sd_age|Subject                          0.07467519                 -0.02634647
#> sigma                                  -0.05632696                  0.01631613
#>                             sd_age|Subject        sigma
#> sd_(Intercept)|Subject         0.074675194 -0.056326961
#> cor_age.(Intercept)|Subject   -0.026346472  0.016316129
#> sd_age|Subject                 0.008375612 -0.004594544
#> sigma                         -0.004594544  0.015888486

if (FALSE) { # \dontrun{
#' # Compare with (parametric) bootstrap results :
get_sdcor <- function(x) {
  as.data.frame(lme4::VarCorr(x), order = "lower.tri")[ , "sdcor"]
}
boo <- bootstrap_mer(fm1, get_sdcor, type = "parametric", nsim = 200L)
# There might be failures in some resamples
cov(boo$t, use = "complete.obs")
} # }
```
