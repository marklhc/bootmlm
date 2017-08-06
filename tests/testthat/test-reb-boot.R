library(nlme)
library(lme4)
library(boot)

m1 <- lmer(Yield ~ 1 | Batch, Dyestuff, REML = FALSE)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
m3 <- lmer(Thickness ~ (1 | Lot) + (1 | Wafer), data = Oxide)

mySumm <- function(.) {
  s <- getME(., "sigma")
  c(beta = getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta")))
}

NSIM <- 19

boo1 <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                      type = "reb", reb_scale = FALSE)
boo1r <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                       type = "reb")
boo1c <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                       type = "reb", reb_scale = TRUE)
boo2 <- bootstrap_mer(m2, mySumm, nsim = NSIM, type = "reb")


# Return Types ------------------------------------------------------------

test_that("random intercept with REB bootstrap", {
  boo <- boo1
  mySumm_m <- mySumm(m1)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
  expect_identical(boo1$t, boo1r$t)
  expect_false(all(boo1$t == boo1c$t))
})

test_that("random slope with REB bootstrap", {
  boo <- boo2
  mySumm_m <- mySumm(m2)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
})

test_that("cross-classified with REB bootstrap", {
  expect_error(bootstrap_mer(m3, mySumm, nsim = NSIM, type = "reb"),
               "REB bootstrap only support one level of clustering")
})


# Bootstrap Confidence Intervals ------------------------------------------

test_that("random intercept with REB bootstrap CI", {
  boo <- boo1
  ci_idx <- sample.int(3, size = 1)
  boo_ci <- boot.ci(boo, index = ci_idx, type = c("norm", "basic", "perc"))
  boo_bca <- boot.ci(boo, index = ci_idx, type = "bca",
                     L = empinf_mer(m1, mySumm, ci_idx))

  expect_output(str(boo_ci), "$ normal", fixed = TRUE)
  expect_output(str(boo_ci), "$ basic", fixed = TRUE)
  expect_output(str(boo_ci), "$ percent", fixed = TRUE)
  expect_error(boot.ci(boo, index = ci_idx, type = "bca"))
  expect_output(str(boo_bca), "$ bca", fixed = TRUE)
  expect_output(str(boo_bca), "num [1, 1:5]", fixed = TRUE)
  expect_warning(boot.ci(boo, index = ci_idx, type = "stud"),
                 "bootstrap variances needed")
})

test_that("random slope with REB bootstrap CI", {
  boo <- boo2
  ci_idx <- sample.int(3, size = 1)
  boo_ci <- boot.ci(boo, index = ci_idx, type = c("norm", "basic", "perc"))
  boo_bca <- boot.ci(boo, index = ci_idx, type = "bca",
                     L = empinf_mer(m2, mySumm, ci_idx))

  expect_output(str(boo_ci), "$ normal", fixed = TRUE)
  expect_output(str(boo_ci), "$ basic", fixed = TRUE)
  expect_output(str(boo_ci), "$ percent", fixed = TRUE)
  expect_error(boot.ci(boo, index = ci_idx, type = "bca"))
  expect_output(str(boo_bca), "$ bca", fixed = TRUE)
  expect_output(str(boo_bca), "num [1, 1:5]", fixed = TRUE)
  expect_warning(boot.ci(boo, index = ci_idx, type = "stud"),
                 "bootstrap variances needed")
})
