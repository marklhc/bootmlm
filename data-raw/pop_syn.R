## Script to generate the `pop_syn` synthetic dataset
## Run this script to reproduce data/pop_syn.rda

library(lme4)
library(haven)
library(dplyr)

# Download and prepare original data
popdata <- haven::read_dta(
  "https://stats.oarc.ucla.edu/stat/stata/examples/mlm_ma_hox/popular.dta"
)

# Fit a random intercept model to extract structural parameters
pop_fit <- lmer(as.numeric(popular) ~ as.numeric(sex) + texp + (1 | school),
                data = popdata, REML = FALSE)

# Estimated parameters (used as the generative model)
# Fixed: intercept = 3.56, sex = 0.84, texp = 0.093
# sigma_u = 0.69, sigma_e = 0.68

# Preserve original school sizes and texp values
school_df <- popdata |>
  group_by(school) |>
  summarise(texp = first(texp), n = n(), .groups = "drop")

n_schools <- nrow(school_df)
sigma_u <- 0.6898
sigma_e <- 0.6780
beta <- c(3.56066, 0.84471, 0.09345)  # intercept, sex, texp

set.seed(2025)
u_j <- rnorm(n_schools, 0, sigma_u)

pop_syn <- purrr::map_dfr(seq_len(n_schools), function(j) {
  nj <- school_df$n[j]
  texp_j <- school_df$texp[j]
  sex_i <- rbinom(nj, 1, 0.5)
  e_ij <- rnorm(nj, 0, sigma_e)
  mu_ij <- beta[1] + beta[2] * sex_i + beta[3] * texp_j + u_j[j] + e_ij
  popular_i <- pmin(pmax(round(mu_ij), 0), 10)
  data.frame(
    pupil = seq_len(nj),
    school = j,
    popular = popular_i,
    sex = sex_i,
    texp = texp_j
  )
})

usethis::use_data(pop_syn, overwrite = TRUE, compress = "bzip2")
