# Score Functions and Case-wise Derivatives

Score Functions and Case-wise Derivatives

## Usage

``` r
scores_mer(x, level = 2)
```

## Arguments

- x:

  A fitted `merMod` object from
  [`lmer`](https://rdrr.io/pkg/lme4/man/lmer.html).

- level:

  If `level = 1`, scores at level-1 are returned; if `level = 2`, which
  is the default, aggregated scores at the cluster- level are returned.
