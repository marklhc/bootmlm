context("Reflating residuals")

library(nlme)
library(lme4)

# traditional representation of BLUP
bstar_trad <- function(x) {
  b <- lme4::getME(x, "b")
  Zt <- lme4::getME(x, "Zt")
  Lambdat <- lme4::getME(x, "Lambdat")
  D <- crossprod(Lambdat)
  X <- model.matrix(x)
  V <- crossprod(lme4::getME(x, "A")) + diag(ncol(Zt))
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

test_that("solve_eigen_sqrt() gives correct result", {
  x <- MASS::mvrnorm(10, mu = c(0, 0, 0),
                     Sigma = matrix(c(1.5, 0.25, -0.2,
                                      0.25, 0.9, 0.02,
                                      -0.2, 0.02, 0.45), nrow = 3))
  M <- var(x)
  expect_equal(tcrossprod(solve_eigen_sqrt(M, M)), M)
})

test_that("same results with get_zeta() and get_zeta_eigen()", {
  M <- Matrix::diag(runif(10))
  r <- rnorm(10)
  expect_equal(get_zeta(r, Matrix::chol(M)), get_zeta_eigen(r, M))
})

test_that("get_reb_resid are correct for RS", {
  Zt <- m1@pp$Zt
  r <- with(m1@resp, y - mu)
  reb_resid <- get_reb_resid(m1, Zt, r, scale = FALSE)
  reb_resid_s <- get_reb_resid(m1, Zt, r, scale = TRUE)
  m1_ml <- reb_resid$ml
  m1_ml_s <- reb_resid_s$ml
  m1_bvc <- tcrossprod(m1_ml_s[[1]]) / ngrps(m1)
  m1_el <- reb_resid$el
  m1_el_s <- reb_resid_s$el

  expect_equal(length(unlist(m1_ml)), length(unlist(m1_ml_s)),
               length(getME(m1, "b")))
  expect_equivalent(m1_bvc, VarCorr(m1)[[1]])
  expect_true(all(diag(tcrossprod(m1_ml_s[[1]])) >
                    diag(tcrossprod(m1_ml[[1]]))))
  expect_equal(length(unlist(m1_el)), length(unlist(m1_el_s)),
               length(resid(m1)))
  expect_true(all(unlist(m1_el_s) / unlist(m1_el) > 1))
})
