library(nlme)
library(lme4)

m1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
m2 <- lmer(Thickness ~ (1 | Lot) + (1 | Wafer), data = Oxide)

test_that(".resid_resample() gives correct output length", {
  boo1 <- .resid_resample(m1, nsim = 3, type = "residual")
  boo2 <- .resid_resample(m2, nsim = 2, type = "residual_cgr")
  boo3 <- .resid_resample(m1, nsim = 3, type = "residual_trans")

  expect_output(str(boo1), "List of 3")
  expect_output(str(boo2), "List of 2")
  expect_output(str(boo3), "num [1:180]", fixed = TRUE)
})

test_that(".case_resample() gives correct output type", {
  boo1 <- .case_resample(m1, nsim = 3, seed = 123, lv1_resample = FALSE)
  boo2 <- .case_resample(m1, nsim = 2, seed = 123, lv1_resample = TRUE)

  expect_output(str(boo1), "List of 3")
  expect_output(str(boo1), "180 obs")
  expect_output(str(boo2), "List of 2")
  expect_output(str(boo2), "3 variables")
  expect_false(identical(boo1[1], boo1[2]))
  expect_false(identical(boo1[1:2], boo2))
})

test_that(".resid_resample() gives distinct resamples", {
  boo1 <- .resid_resample(m1, nsim = 2, seed = 122, type = "residual")
  boo2 <- .resid_resample(m1, nsim = 3, seed = 122, type = "residual")
  boo3 <- .resid_resample(m1, nsim = 2, seed = 123, type = "residual")
  boo4 <- .resid_resample(m1, seed = 123, type = "residual_trans",
                          corrected = FALSE)
  boo5 <- .resid_resample(m1, seed = 123, type = "residual_trans",
                          corrected = TRUE)

  expect_false(identical(boo1[1], boo1[2]))
  expect_identical(boo1, boo2[1:2])
  expect_false(identical(boo1, boo3))
  expect_false(identical(boo4, boo5))
})
