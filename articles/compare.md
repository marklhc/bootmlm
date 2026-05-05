# Comparing Different Multilevel Bootstrapping Methods

## Data

Here I use a synthetic dataset (`pop_syn`) bundled with the package,
which mimics the structure of the example data in Hox (2010, p. 17). The
synthetic data was generated from parameters estimated from the original
*popular* dataset.

``` r

library(dplyr)
library(lme4)
library(boot)
library(bootmlm)
```

For the interest of time, I will select only 30 schools and a few
students in each school:

``` r

set.seed(85957)
pop_sub <- pop_syn |> filter(school %in% sample(unique(school), 30)) |>
  group_by(school) |> sample_frac(size = .25) |> ungroup()
```

To illustrate doing bootstrapping to obtain standard error (*SE*) and
confidence interval (CI) for the intraclass correlation (ICC) of the
variable `popular`, we first fit an intercept only model:

``` r

m0 <- lmer(popular ~ (1 | school), data = pop_sub)
summary(m0)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: popular ~ (1 | school)
    ##    Data: pop_sub
    ## 
    ## REML criterion at convergence: 441.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.84638 -0.60538 -0.03404  0.67043  2.42008 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  school   (Intercept) 0.9692   0.9845  
    ##  Residual             0.7110   0.8432  
    ## Number of obs: 152, groups:  school, 30
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept)   5.2634     0.1925   27.34

The intraclass correlation is defined as
``` math
\rho = \frac{\tau}{\tau + \sigma^2} = \frac{1}{1 + \sigma^2 / \tau} 
  = \frac{1}{1 + \theta^{-2}},
```
where $`\theta = \sqrt{\tau} / \sigma`$ is the relative cholesky factor
for the random intercept term used in `lme4`. Therefore, we can estimate
the ICC as:

``` r

(icc0 <- 1 / (1 + getME(m0, "theta")^(-2)))
```

    ## school.(Intercept) 
    ##          0.5768188

So the ICC is quite large for this data set. However, it is important to
also quantify the uncertainty of a point estimate. Although there are
analytic methods to obtain *SE* and CI for ICC, a reliable alternative
is to do bootstrapping. We first define the function for computing the
test statistic:

``` r

icc <- function(x) 1 / (1 + x@theta^(-2))
icc(m0)
```

    ## [1] 0.5768188

Note that I used the `@` operator for faster extraction. With the
`bootmlm` package we can perform various bootstrap methods using the
[`bootstrap_mer()`](https://marklhc.github.io/bootmlm/reference/bootstrap_mer.md)
function.

## Parametric Bootstrap

We can run parametric bootstrap, which essentially call the
[`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html) function.
It’s usually recommended to have a large number of bootstrap samples
($`R`$), especially for CIs with higher confidence levels. For
illustrative purpose I will use $`R = 999`$, but in general 1,999 or
more is recommended

``` r

system.time(
  boo_par <- bootstrap_mer(m0, icc, nsim = 999L, type = "parametric")
)
boo_par
```

    ## 
    ## PARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## lme4::bootMer(x = x, FUN = FUN, nsim = nsim, seed = seed, use.u = FALSE, 
    ##     type = "parametric", verbose = FALSE)
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original      bias    std. error
    ## t1* 0.5768188 -0.01019232  0.08598194

As you can see, the *SE* for the ICC is estimated to be 0.0859819 with
parametric bootstrap.

### Confidence interval

With parametric bootstrap there are three ways to construct confidence
intervals via the
[`boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) function from
the `boot` package: normal, basic, and percentile. We can use the
following function:

``` r

boo_par_ci <- boot::boot.ci(boo_par, type = c("norm", "basic", "perc"), 
                            index = 1L)
boo_par_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_par, type = c("norm", "basic", "perc"), 
    ##     index = 1)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic              Percentile     
    ## 95%   ( 0.4185,  0.7555 )   ( 0.4369,  0.7813 )   ( 0.3723,  0.7167 )  
    ## Calculations and Intervals on Original Scale

## Residual Bootstraps

Whereas parametric bootstrap resamples from independent normal
distributions, residual bootstrap samples the residuals. Therefore,
residual bootstrap is expected to be more robust to non-normality.
`bootmlm` implements three methods for residual bootstrap:
differentially reflated residual bootstrap,
Carpenter-Goldstein-Rashbash’s residual bootstrap (CGR; Carpenter et
al., 2003), and transformational residual bootstrap by van der Leeden,
Meijer, and Busin (2008). They are all motivated by the fact that the
residuals, generally empirical bayes estimates (denoted as $`\tilde u`$
and $`\tilde e`$), are shrinkage estimates and have sampling
variabilities much smaller than the population random effects, $`u`$ and
$`e`$.

The first residual bootstrap rescale $`\tilde u`$ and $`\tilde e`$ so
that their sampling variabilities match those of $`u`$ and $`e`$ as
implied by the model estimates. This can be obtained by

``` r

system.time(
  boo_res <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual")
)
boo_res
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "residual")
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original       bias    std. error
    ## t1* 0.5768188 -0.008674732  0.07221562

As you can see, the *SE* for the ICC is estimated to be 0.0722156 with
residual bootstrap.

The second method, CGR, rescale the sample covariance matrix of the
*realized values* of the residuals to match the model-implied variance
components. This can be obtained by

``` r

system.time(
  boo_cgr <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual_cgr")
)
boo_cgr
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "residual_cgr")
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original       bias    std. error
    ## t1* 0.5768188 -0.007783322  0.07040576

The *SE* is estimated to be 0.0704058 with CGR bootstrap.

The third method first transforms the OLS residuals,
$`\hat r_{ij} = y_{ij} -
\boldsymbol{x}_{ij} \hat{\boldsymbol{\beta}}`$, by the inverse of
cholesky factor, $`\boldsymbol{L}`$, of the model-implied covariance
matrix of $`\boldsymbol{y}`$, $`\hat{\boldsymbol{V}}`$, so that
theoretically
$`\boldsymbol{L}^{-1} (\boldsymbol{y} - \boldsymbol{X \beta})`$ should
be independent and identically distributed. However, as the true
sampling variance of $`\hat r_{ij}`$ is not $`\boldsymbol{V}`$, I also
provide the option `corrected_trans = TRUE` to do the transformation
using the theoretically sampling variability of $`\hat r_{ij}`$.

``` r

# Transformation according to V
system.time(
  boo_tra <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual_trans")
)
boo_tra

# Transformation according to the sampling variance of r
system.time(
  boo_trac <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual_trans", 
                            corrected_trans = TRUE)
)
boo_trac
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "residual_trans")
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original       bias    std. error
    ## t1* 0.5768188 -0.009255973  0.07836078

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "residual_trans", 
    ##     corrected_trans = TRUE)
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original      bias    std. error
    ## t1* 0.5768188 -0.01245717  0.08085008

The *SE* is estimated to be 0.0783608 and 0.0808501 with and without
corrections with the transformational residual bootstrap.

### Confidence interval

With residual bootstrap methods there are four ways to construct
confidence intervals via the
[`boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) function from
the `boot` package, with the addition of the bias-corrected and
accelarted bootstrap (BCa). We can use the following function:

``` r

# First need to compute the influence values
inf_val <- empinf_mer(m0, icc, index = 1)

# Residual bootstrap
boo_res_ci <- boot::boot.ci(boo_res, type = c("norm", "basic", "perc", "bca"), 
                            index = 1L, L = inf_val)
boo_res_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_res, type = c("norm", "basic", "perc", 
    ##     "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.4440,  0.7270 )   ( 0.4583,  0.7456 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.4080,  0.6953 )   ( 0.4465,  0.7146 )  
    ## Calculations and Intervals on Original Scale
    ## Some BCa intervals may be unstable

``` r

# CGR
boo_cgr_ci <- boot::boot.ci(boo_cgr, type = c("norm", "basic", "perc", "bca"), 
                            index = 1L, L = inf_val)
boo_cgr_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_cgr, type = c("norm", "basic", "perc", 
    ##     "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.4466,  0.7226 )   ( 0.4538,  0.7281 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.4255,  0.6998 )   ( 0.4472,  0.7244 )  
    ## Calculations and Intervals on Original Scale
    ## Some BCa intervals may be unstable

``` r

# Transformational (no correction)
boo_tra_ci <- boot::boot.ci(boo_tra, type = c("norm", "basic", "perc", "bca"), 
                            index = 1L, L = inf_val)
boo_tra_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_tra, type = c("norm", "basic", "perc", 
    ##     "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.4325,  0.7397 )   ( 0.4570,  0.7664 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.3872,  0.6966 )   ( 0.4241,  0.7081 )  
    ## Calculations and Intervals on Original Scale

``` r

# Transformational (with correction)
boo_trac_ci <- boot::boot.ci(boo_trac, type = c("norm", "basic", "perc", "bca"), 
                             index = 1L, L = inf_val)
boo_trac_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_trac, type = c("norm", "basic", 
    ##     "perc", "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.4308,  0.7477 )   ( 0.4498,  0.7624 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.3913,  0.7038 )   ( 0.4229,  0.7331 )  
    ## Calculations and Intervals on Original Scale
    ## Some BCa intervals may be unstable

## Random Effect Block Bootstrap

``` r

system.time(
  boo_reb <- bootstrap_mer(m0, icc, nsim = 999L, type = "reb")
)
boo_reb

system.time(
  boo_rebs <- bootstrap_mer(m0, icc, nsim = 999L, type = "reb", 
                            reb_scale = TRUE)
)
boo_rebs
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "reb")
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original     bias    std. error
    ## t1* 0.5768188 0.06957508  0.06176538

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "reb", reb_scale = TRUE)
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original      bias    std. error
    ## t1* 0.5768188 -0.01103587  0.07066787

### Confidence interval

``` r

# Only sampling clusters
boo_reb_ci <- boot::boot.ci(boo_reb, type = c("norm", "basic", "perc", "bca"), 
                            index = 1L, L = inf_val)
```

    ## Warning in norm.inter(t, adj.alpha): extreme order statistics used as endpoints

``` r

boo_reb_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_reb, type = c("norm", "basic", "perc", 
    ##     "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.3862,  0.6283 )   ( 0.4026,  0.6408 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.5129,  0.7510 )   ( 0.3551,  0.6389 )  
    ## Calculations and Intervals on Original Scale
    ## Warning : BCa Intervals used Extreme Quantiles
    ## Some BCa intervals may be unstable

``` r

# Transformational (with correction)
boo_rebs_ci <- boot::boot.ci(boo_rebs, type = c("norm", "basic", "perc", "bca"), 
                             index = 1L, L = inf_val)
boo_rebs_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_rebs, type = c("norm", "basic", 
    ##     "perc", "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.4493,  0.7264 )   ( 0.4642,  0.7419 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.4117,  0.6895 )   ( 0.4521,  0.7157 )  
    ## Calculations and Intervals on Original Scale
    ## Some BCa intervals may be unstable

## Case Bootstrap

With case bootstrap, the observed *cases* are sampled with replacement.
However, because of the multilevel structure, we need to resample the
clusters. Optionally, we can then resample the cases within each cluster
(using the `lv1_resample = TRUE` argument). Unlike the parametric and
residual bootstrap methods, currently `bootmlm` only support the case
bootstrap with two levels.

``` r

system.time(
  boo_cas <- bootstrap_mer(m0, icc, nsim = 999L, type = "case")
)
boo_cas

system.time(
  boo_cas1 <- bootstrap_mer(m0, icc, nsim = 999L, type = "case", 
                            lv1_resample = TRUE)
)
boo_cas1
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "case")
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original      bias    std. error
    ## t1* 0.5768188 -0.01577482  0.05671946

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## bootstrap_mer(x = m0, FUN = icc, nsim = 999, type = "case", lv1_resample = TRUE)
    ## 
    ## 
    ## Bootstrap Statistics :
    ##      original     bias    std. error
    ## t1* 0.5768188 0.06758788  0.06112113

The *SE* for the ICC is estimated to be 0.0567195 (only sampling
clusters) and 0.0611211 (sampling also cases) with case bootstrap.

### Confidence interval

With case bootstrap the supported CIs are: normal, basic, and
percentile, and BCa. We can use the following function:

``` r

# Only sampling clusters
boo_cas_ci <- boot::boot.ci(boo_cas, type = c("norm", "basic", "perc", "bca"), 
                            index = 1L, L = inf_val)
boo_cas_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_cas, type = c("norm", "basic", "perc", 
    ##     "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.4814,  0.7038 )   ( 0.4953,  0.7223 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.4314,  0.6583 )   ( 0.4807,  0.6838 )  
    ## Calculations and Intervals on Original Scale
    ## Some BCa intervals may be unstable

``` r

# Transformational (with correction)
boo_cas1_ci <- boot::boot.ci(boo_cas1, type = c("norm", "basic", "perc", "bca"), 
                             index = 1L, L = inf_val)
```

    ## Warning in norm.inter(t, adj.alpha): extreme order statistics used as endpoints

``` r

boo_cas1_ci
```

    ## BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
    ## Based on 999 bootstrap replicates
    ## 
    ## CALL : 
    ## boot::boot.ci(boot.out = boo_cas1, type = c("norm", "basic", 
    ##     "perc", "bca"), index = 1, L = inf_val)
    ## 
    ## Intervals : 
    ## Level      Normal              Basic         
    ## 95%   ( 0.3894,  0.6290 )   ( 0.4008,  0.6435 )  
    ## 
    ## Level     Percentile            BCa          
    ## 95%   ( 0.5101,  0.7529 )   ( 0.4189,  0.6407 )  
    ## Calculations and Intervals on Original Scale
    ## Warning : BCa Intervals used Extreme Quantiles
    ## Some BCa intervals may be unstable

## Summary

``` r

boo_names <- c("parametric", "residual", "cgr", "trans", 
               "trans (cor)", "REB", "REB (scaled)", 
               "case (cluster)", "case (c + i)")
boo_lst <- list(boo_par, boo_res, boo_cgr, boo_tra, boo_trac, 
                boo_reb, boo_rebs, boo_cas, boo_cas1)
boo_ci_lst <- list(boo_par_ci, boo_res_ci, boo_cgr_ci, boo_tra_ci, 
                   boo_trac_ci, boo_reb_ci, boo_rebs_ci, boo_cas_ci, 
                   boo_cas1_ci)

get_ci <- function(boo_ci, type) {
  paste0("(", paste(comma(tail(boo_ci[[type]][1, ], 2L)), collapse = ", "), ")")
}

tab <- tibble(boot_type = boo_names, boo = boo_lst, boo_ci = boo_ci_lst) |>
  mutate(sd = sapply(boo, \(x) comma(sd(x$t))), 
         normal = sapply(boo_ci, \(x) get_ci(x, "normal")), 
         basic = sapply(boo_ci, \(x) get_ci(x, "basic")), 
         percentile = sapply(boo_ci, \(x) get_ci(x, "percent")), 
         bca = sapply(boo_ci, \(x) get_ci(x, "bca"))) |>
  select(-boo, -boo_ci)

knitr::kable(tab)
```

| boot_type      | sd    | normal       | basic        | percentile   | bca          |
|:---------------|:------|:-------------|:-------------|:-------------|:-------------|
| parametric     | 0.086 | (0.42, 0.76) | (0.44, 0.78) | (0.37, 0.72) | (NULL)       |
| residual       | 0.072 | (0.44, 0.73) | (0.46, 0.75) | (0.41, 0.70) | (0.45, 0.71) |
| cgr            | 0.07  | (0.45, 0.72) | (0.45, 0.73) | (0.43, 0.70) | (0.45, 0.72) |
| trans          | 0.078 | (0.43, 0.74) | (0.46, 0.77) | (0.39, 0.70) | (0.42, 0.71) |
| trans (cor)    | 0.081 | (0.43, 0.75) | (0.45, 0.76) | (0.39, 0.70) | (0.42, 0.73) |
| REB            | 0.062 | (0.39, 0.63) | (0.40, 0.64) | (0.51, 0.75) | (0.36, 0.64) |
| REB (scaled)   | 0.071 | (0.45, 0.73) | (0.46, 0.74) | (0.41, 0.69) | (0.45, 0.72) |
| case (cluster) | 0.057 | (0.48, 0.70) | (0.50, 0.72) | (0.43, 0.66) | (0.48, 0.68) |
| case (c + i)   | 0.061 | (0.39, 0.63) | (0.40, 0.64) | (0.51, 0.75) | (0.42, 0.64) |
