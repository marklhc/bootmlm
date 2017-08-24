library(lme4)

m1 <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)

test_that("empinf_mer() and empinf_merm() give same results", {
  L1 <- empinf_mer(m1, sigma)
  L2 <- empinf_merm(m1, function(x) c(x@theta, sigma(x)))

  expect_length(L1, lme4::ngrps(m1))
  expect_equal(ncol(L2), 2)
  expect_identical(L1, L2[ , 2])
})
