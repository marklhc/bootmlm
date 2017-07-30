# traditional representation of BLUP
library(nlme)
library(lme4)
bstar_trad <- function(x) {
  y <- getME(x, "y")
  b <- getME(x, "b")
  Zt <- getME(x, "Zt")
  Lambdat <- getME(x, "Lambdat")
  D <- crossprod(Lambdat)
  X <- model.matrix(x)
  V <- crossprod(getME(x, "A")) + diag(ncol(Zt))
  VinvX <- Matrix::solve(V, X)
  Vb <- D %*% Zt %*% Matrix::solve(V, crossprod(Zt, D)) -
    D %*% Zt %*% VinvX %*% Matrix::solve(crossprod(X, VinvX),
                                         crossprod(VinvX, crossprod(Zt, D)))
  Vb[D == 0] <- 0
  bstar <- crossprod(Lambdat, Matrix::solve(t(base::chol(Vb)), b))
  bstar
}

m1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
m2 <- lmer(Thickness ~ (1 | Lot) + (1 | Wafer), data = Oxide)

test_that("get_reflate_b() and get_reflate_e() are correct for RS", {
  m1_e <- resid(m1)
  m1_estar <- get_reflate_e(m1)
  m1_bstar <- get_reflate_b(m1)

  expect_length(unlist(m1_bstar), length(getME(m1, "b")))
  expect_equal(unlist(m1_bstar), as.vector(bstar_trad(m1)))
  expect_length(m1_estar, length(m1_e))
  expect_true(all(m1_estar / m1_e > 1))
})

test_that("get_reflate_b() and get_reflate_e() are correct for CCREM", {
  m2_e <- resid(m2)
  m2_estar <- get_reflate_e(m2)
  m2_bstar <- get_reflate_b(m2)

  expect_length(unlist(m2_bstar), length(getME(m2, "b")))
  expect_equal(unlist(m2_bstar), as.vector(bstar_trad(m2)))
  expect_length(m2_estar, length(m2_e))
  expect_true(all(m2_estar / m2_e > 1))
})

test_that("get_reflate_b_cgr() gives reflated residuals", {
  m1_bstar <- get_reflate_b_cgr(m1)
  m1_bvc <- tcrossprod(m1_bstar[[1]]) / ngrps(m1)

  m2_bstar <- get_reflate_b_cgr(m2)
  m2_bvc <- lapply(seq_along(m2_bstar),
                   function(i) tcrossprod(m2_bstar[[i]]) / ngrps(m2)[[i]])

  expect_equivalent(m1_bvc, VarCorr(m1)[[1]])
  expect_equivalent(m2_bvc, VarCorr(m2))
})
