# bootmlm

The `bootmlm` package does bootstrap resampling for multilevel models. 
Currently only models fitted with `lme4::lmer()` is supported. It's still in 
developmental stage and is not yet on CRAN. However, you can install the package on GitHub:

```r
remotes::install_github("marklhc/bootmlm")
```

See [this paper](https://quantscience.rbind.io/files/Lai_2020_mbr_bootstrap_mlm_effectsize_am.pdf) for a performance comparison of different bootstrapped confidence
intervals for multilevel effect size estimations:

Lai, M. H. C. (2020). Bootstrap confidence interval for multilevel standardized effect size. Multivariate Behavioral Research. Advance online publication. https://doi.org/10.1080/00273171.2020.1746902

## Example

Here is an example to get the bootstrap distributions of the fixed effects and
the level-1 error SD:

```r
library(lme4)
fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
mySumm <- function(x) {
  c(getME(x, "beta"), sigma(x))
}
# Covariance preserving residual bootstrap
library(bootmlm)
boo01 <- bootstrap_mer(fm01ML, mySumm, type = "residual", nsim = 100)
# Plot bootstrap distribution of fixed effect
library(boot)
plot(boo01, index = 1)
# Get confidence interval
boot.ci(boo01, index = 2, type = c("norm", "basic", "perc"))
# BCa using influence values computed from `empinf_mer`
boot.ci(boo01, index = 2, type = "bca", L = empinf_mer(fm01ML, mySumm, 2))
```
