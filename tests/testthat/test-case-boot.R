context("Case bootstrap")

library(nlme)
library(lme4)
library(boot)

m1 <- lmer(Yield ~ 1 | Batch, Dyestuff, REML = FALSE)
m2 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
m3 <- lmer(Thickness ~ (1 | Lot) + (1 | Wafer), data = Oxide)
m4 <- lmer(angle ~ recipe * temperature + (1 | recipe:replicate), cake,
           REML= FALSE)

mySumm <- function(.) {
  s <- lme4::getME(., "sigma")
  c(beta = lme4::getME(., "beta"), sigma = s,
    sig01 = unname(s * lme4::getME(., "theta")))
}

NSIM <- 19

boo1 <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                      type = "case", lv1_resample = FALSE)
boo1r <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                       type = "case")
boo1c <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                       type = "case", lv1_resample = TRUE)
boo2 <- bootstrap_mer(m2, mySumm, nsim = NSIM, type = "case")
boo4 <- bootstrap_mer(m4, mySumm, nsim = NSIM, type = "case")


# Return Types ------------------------------------------------------------

test_that("random intercept with case bootstrap", {
  boo <- boo1
  mySumm_m <- mySumm(m1)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
  expect_identical(boo1$t, boo1r$t)
  expect_false(all(boo1$t == boo1c$t))
})

test_that("random slope with case bootstrap", {
  boo <- boo2
  mySumm_m <- mySumm(m2)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
})

test_that("cross-classified with case bootstrap", {
  expect_error(bootstrap_mer(m3, mySumm, nsim = NSIM, type = "case"),
               "case bootstrap only support one level of clustering")
})

test_that("interaction with case bootstrap", {
  boo <- boo4
  mySumm_m <- mySumm(m4)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
})


# Bootstrap Confidence Intervals ------------------------------------------

test_that("random intercept with case bootstrap CI", {
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

test_that("random slope with case bootstrap CI", {
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

test_that("interaction with case bootstrap CI", {
  boo <- boo4
  ci_idx <- sample.int(5, size = 1)
  boo_ci <- boot.ci(boo, index = ci_idx, type = c("norm", "basic", "perc"))
  boo_bca <- boot.ci(boo, index = ci_idx, type = "bca",
                     L = empinf_mer(m4, mySumm, ci_idx))

  expect_output(str(boo_ci), "$ normal", fixed = TRUE)
  expect_output(str(boo_ci), "$ basic", fixed = TRUE)
  expect_output(str(boo_ci), "$ percent", fixed = TRUE)
  expect_error(boot.ci(boo, index = ci_idx, type = "bca"))
  expect_output(str(boo_bca), "$ bca", fixed = TRUE)
  expect_output(str(boo_bca), "num [1, 1:5]", fixed = TRUE)
  expect_warning(boot.ci(boo, index = ci_idx, type = "stud"),
                 "bootstrap variances needed")
})
