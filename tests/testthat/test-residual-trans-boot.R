context("transformational residual bootstrap")

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
                      type = "residual_trans")
boo1r <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                       type = "residual_trans")
boo1c <- bootstrap_mer(m1, mySumm, nsim = NSIM, seed = 101,
                       corrected_trans = TRUE,
                       type = "residual_trans")
boo2 <- bootstrap_mer(m2, mySumm, nsim = NSIM, type = "residual_trans")
boo3 <- bootstrap_mer(m3, mySumm, nsim = NSIM, type = "residual_trans",
                      corrected_trans = TRUE)
boo4 <- bootstrap_mer(m4, mySumm, nsim = NSIM, type = "residual_trans")


# Return Types ------------------------------------------------------------

test_that("random intercept with transformed residual bootstrap", {
  boo <- boo1
  mySumm_m <- mySumm(m1)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
  expect_identical(boo1$t, boo1r$t)
  expect_false(all(boo1$t == boo1c$t))
})

test_that("random slope with transformed residual bootstrap", {
  boo <- boo2
  mySumm_m <- mySumm(m2)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
})

test_that("cross-classified with transformed residual bootstrap", {
  boo <- boo3
  mySumm_m <- mySumm(m3)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
})

test_that("interaction with transformed residual bootstrap", {
  boo <- boo4
  mySumm_m <- mySumm(m4)

  expect_s3_class(boo, "boot")
  expect_equal(boo$t0, mySumm_m)
  expect_equal(nrow(boo$t), boo$R, NSIM)
  expect_equal(ncol(boo$t), length(mySumm_m))
})


# Bootstrap Confidence Intervals ------------------------------------------

test_that("random intercept with transformed residual bootstrap CI", {
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

test_that("random slope with transformed residual bootstrap CI", {
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

test_that("cross-classified with transformed residual bootstrap CI", {
  boo <- boo3
  ci_idx <- sample.int(3, size = 1)
  boo_ci <- boot.ci(boo, index = ci_idx, type = c("norm", "basic", "perc"))

  expect_output(str(boo_ci), "$ normal", fixed = TRUE)
  expect_output(str(boo_ci), "$ basic", fixed = TRUE)
  expect_output(str(boo_ci), "$ percent", fixed = TRUE)
  expect_error(boot.ci(boo, index = ci_idx, type = "bca"))
  expect_error(boot.ci(boo, index = ci_idx, type = "bca",
                       L = empinf_mer(m3, mySumm, ci_idx)))
  expect_warning(boot.ci(boo, index = ci_idx, type = "stud"),
                 "bootstrap variances needed")
})

test_that("interaction with transformed residual bootstrap CI", {
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
