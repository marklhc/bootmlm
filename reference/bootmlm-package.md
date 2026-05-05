# bootmlm: A package for bootstrap resampling with multilevel data.

Currently, the [`bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html)
function in the lme4 package only implements the parametric bootstrap
and a limited version of semiparametric (or residual) bootstrap. The
bootmlm package provides the function `bootstrap_mer`, which performs
various parametric, residual, and case bootstrap resampling for fitted
model objects with the lme4 package.

## Limitations

- Currently only support multilevel models (a.k.a. linear mixed-effects
  models) fitted by [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)
  with the lme4 package. Support for categorical outcome (i.e.,
  generalized linear mixed-effects models fitted by
  [`glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)) and for models
  fitted with the nlme package will be added in the future.

- Random effect block bootstrap (`type = 'reb'`) and case bootstrap
  (`type = 'case'`) only support two-level models.

- Bias-corrected and accelerated bootstrap (using
  [`empinf_mer()`](https://marklhc.github.io/bootmlm/reference/empinf_mer.md))
  only supports two-level models.

## See also

Useful links:

- <https://github.com/marklhc/bootmlm>

- Report bugs at <https://github.com/marklhc/bootmlm/issues>

## Author

**Maintainer**: Mark Lai <marklhc@gmail.com>

Authors:

- Mark Lai <marklhc@gmail.com>
