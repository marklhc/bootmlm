# data-raw/vignette_data_derived.R
# Run this locally, NEVER on CRAN.

library(dplyr)
library(lme4)
library(boot)
library(bootmlm)
library(msm)

# 1. Setup the exact data from the vignette
set.seed(85957)
pop_sub <- pop_syn |> filter(school %in% sample(unique(school), 20)) |>
  group_by(school) |> sample_frac(size = .25) |> ungroup()

m0 <- lmer(popular ~ (1 | school), data = pop_sub)

# 2. Define the extraction function
icc <- function(x) {
  th_est <- x@theta
  est <- 1 / (1 + th_est^(-2))
  th_var <- vcov_theta(x)
  var <- msm::deltamethod(g = ~ 1 / (1 + x1^(-2)), 
                            mean = th_est, cov = th_var, ses = FALSE)
  c(est, var)
}

# 3. Run the heavy computations
boo <- bootstrap_mer(m0, icc, 999L, type = "case")
inf_val <- empinf_mer(m0, icc, index = 1)

# 4. Save to the vignettes folder
save(boo, inf_val, file = "vignettes/precomputed_derived.rda", compress = "xz")