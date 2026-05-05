# data-raw/vignette_computations.R
# Run this locally, NEVER on CRAN.

library(dplyr)
library(lme4)
library(boot)
library(bootmlm)

# 1. Setup the data exactly as in the vignette
set.seed(85957)
pop_sub <- pop_syn |> filter(school %in% sample(unique(school), 30)) |>
  group_by(school) |> sample_frac(size = .25) |> ungroup()

m0 <- lmer(popular ~ (1 | school), data = pop_sub)
icc <- function(x) 1 / (1 + x@theta^(-2))

# 2. Run all the heavy bootstraps
boo_par  <- bootstrap_mer(m0, icc, nsim = 999L, type = "parametric")
boo_res  <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual")
boo_cgr  <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual_cgr")
boo_tra  <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual_trans")
boo_trac <- bootstrap_mer(m0, icc, nsim = 999L, type = "residual_trans", corrected_trans = TRUE)
boo_reb  <- bootstrap_mer(m0, icc, nsim = 999L, type = "reb")
boo_rebs <- bootstrap_mer(m0, icc, nsim = 999L, type = "reb", reb_scale = TRUE)
boo_cas  <- bootstrap_mer(m0, icc, nsim = 999L, type = "case")
boo_cas1 <- bootstrap_mer(m0, icc, nsim = 999L, type = "case", lv1_resample = TRUE)

# 3. Save them to the vignettes folder so the .Rmd can find them during the CRAN build
save(boo_par, boo_res, boo_cgr, boo_tra, boo_trac, boo_reb, boo_rebs, boo_cas, boo_cas1, 
     file = "vignettes/precomputed_bootstraps.rda", compress = "xz")