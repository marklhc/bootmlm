var_blup <- function(u, L, RX, A, X) {
  q1 <- Matrix::solve(L, A)
  q2 <- X %*% Matrix::chol2inv(RX) %*% t(X)
  Q <- q1 + q1 %*% q2 %*% crossprod(A, q1) - q1 %*% q2
  # tcrossprod(Q) + crossprod(tcrossprod(A, Q))
  Matrix::rowSums(Q^2) + Matrix::colSums(tcrossprod(Q, A)^2)
}

get_reflate_b <- function(x) {
  # Extract required quantities from the S4 object
  PR <- x@pp
  rsp <- x@resp
  L <- PR$L()
  RX <- PR$RX()
  Lambdat <- PR$Lambdat
  A <- Lambdat %*% PR$Zt
  X <- model.matrix(x)
  u <- x@u

  # Hat matrix for u
  if (all(u == 0)) {  # may break down when variance is 0 for just one component
    ustar <- u
  } else {
    Vu <- var_blup(u, L, RX, A, X)
    ustar <- u / sqrt(Vu)
  }
  crossprod(Lambdat, ustar)
}

get_reflate_e <- function(x) {
  rsp <- x@resp
  hat_e <- hatvalues(x, fullHatMatrix = FALSE)
  (rsp$y - rsp$mu) / sqrt(1 - hat_e)
}
