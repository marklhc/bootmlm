context("Asymptotic covariances for variance components")

library(nlme)
library(lme4)

m1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
m1_lme <- lme(Reaction ~ Days, random = ~ Days | Subject, data = sleepstudy)
m2 <- lmer(Thickness ~ (1 | Lot) + (1 | Wafer), data = Oxide)

vc_fun <- function(x) {
  vdd <- as.data.frame(lme4::VarCorr(x), order = "lower.tri")
  vdd[ , "sdcor"]
}

test_that("vcov_vc() gives reasonable estimates as parametric bootstrap", {
  sd_vc1 <- sqrt(diag(vcov_vc(m1)))
  sd_vc1_lme <- c(with(intervals(m1_lme)$reStruct[[1]],
                         (upper - lower) / 2 / qnorm(.975)),
                    (intervals(m1_lme)$sigma["upper"] -
                       intervals(m1_lme)$sigma["lower"]) / 2 / qnorm(.975))
  sd_vc1_lme <- sd_vc1_lme[c(1, 3, 2, 4)]
  expect_true(all(abs(log(sd_vc1 / sd_vc1_lme)) < 0.1))
})
