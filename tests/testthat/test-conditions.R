library(lme4)
library(boot)

m1 <- lmer(Yield ~ 1 | Batch, Dyestuff, REML = FALSE)
m2 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
            data = cbpp, family = binomial)

mySumm <- function(.) {
  s <- lme4::getME(., "sigma")
  c(beta = lme4::getME(., "beta"), sigma = s,
    sig01 = unname(s * lme4::getME(., "theta")))
}

NSIM <- 19

test_that("bootstrap_mer(..., type = 'parametric') same as BootMer()", {
  boo1 <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                        type = "parametric")
  boo2 <- bootMer(m1, mySumm, nsim = NSIM, seed = 101)
  expect_identical(boo1$t, boo2$t)
})

test_that("incorrect input returns error", {
  boo <- bootstrap_mer(m1, FUN = "mean", nsim = NSIM, type = "residual")
  fun_lst <- function(x) list(sigma = sigma(x))

  expect_true(all(is.na(boo$t)))
  expect_error(bootstrap_mer(m1, FUN = mySumm, nsim = -1L, type = "case"))
  expect_error(
    bootstrap_mer(m2, FUN = mySumm, nsim = NSIM, type = "residual_cgr"),
    "only linear mixed model of class")
  expect_error(
    bootstrap_mer(m1, FUN = fun_lst, nsim = NSIM, type = "residual_trans"),
    "functions that return numeric vectors")
  expect_error(
    bootstrap_mer(m1, FUN = fun_lst, nsim = NSIM, type = "resid"),
    "should be one of")
})

test_that("`corrected_trans` only affect `type = 'residual_trans'", {
  identical_rerun <- function(x, type, FUN = mySumm, nsim = 5L, seed = 123) {
    # Whether output is the same with `corrected_trans` as TRUE or FALSE
    boo1 <- bootstrap_mer(x, FUN = FUN, nsim = nsim, type = type, seed = seed,
                          corrected_trans = FALSE)
    boo2 <- bootstrap_mer(x, FUN = FUN, nsim = nsim, type = type, seed = seed,
                          corrected_trans = TRUE)
    boo1 <- boo1[names(boo1) != "call"]
    boo2 <- boo2[names(boo2) != "call"]
    identical(boo1, boo2)
  }

  expect_true(identical_rerun(m1, "parametric"))
  expect_true(identical_rerun(m1, "residual"))
  expect_true(identical_rerun(m1, "residual_cgr"))
  expect_false(identical_rerun(m1, "residual_trans"))
  expect_true(identical_rerun(m1, "reb"))
  expect_true(identical_rerun(m1, "case"))
})

test_that("`reb_scale` only affect `type = 'reb'", {
  identical_rerun <- function(x, type, FUN = mySumm, nsim = 5L, seed = 124) {
    # Whether output is the same with `corrected_trans` as TRUE or FALSE
    boo1 <- bootstrap_mer(x, FUN = FUN, nsim = nsim, type = type, seed = seed,
                          reb_scale = FALSE)
    boo2 <- bootstrap_mer(x, FUN = FUN, nsim = nsim, type = type, seed = seed,
                          reb_scale = TRUE)
    boo1 <- boo1[names(boo1) != "call"]
    boo2 <- boo2[names(boo2) != "call"]
    identical(boo1, boo2)
  }

  expect_true(identical_rerun(m1, "parametric"))
  expect_true(identical_rerun(m1, "residual"))
  expect_true(identical_rerun(m1, "residual_cgr"))
  expect_true(identical_rerun(m1, "residual_trans"))
  expect_false(identical_rerun(m1, "reb"))
  expect_true(identical_rerun(m1, "case"))
})


test_that("`lv1_resample` only affect `type = 'case'", {
  identical_rerun <- function(x, type, FUN = mySumm, nsim = 5L, seed = 124) {
    # Whether output is the same with `corrected_trans` as TRUE or FALSE
    boo1 <- bootstrap_mer(x, FUN = FUN, nsim = nsim, type = type, seed = seed,
                          lv1_resample = FALSE)
    boo2 <- bootstrap_mer(x, FUN = FUN, nsim = nsim, type = type, seed = seed,
                          lv1_resample = TRUE)
    boo1 <- boo1[names(boo1) != "call"]
    boo2 <- boo2[names(boo2) != "call"]
    identical(boo1, boo2)
  }

  expect_true(identical_rerun(m1, "parametric"))
  expect_true(identical_rerun(m1, "residual"))
  expect_true(identical_rerun(m1, "residual_cgr"))
  expect_true(identical_rerun(m1, "residual_trans"))
  expect_true(identical_rerun(m1, "reb"))
  expect_false(identical_rerun(m1, "case"))
})
