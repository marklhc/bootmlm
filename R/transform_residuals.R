var_blup <- function(L, RX, A, Lambdat, X) {
  q1 <- Matrix::solve(L, A)
  q2 <- X %*% Matrix::chol2inv(RX) %*% t(X)
  Q <- q1 + q1 %*% q2 %*% crossprod(A, q1) - q1 %*% q2
  LambdaQ <- crossprod(Lambdat, Q)
  tcrossprod(LambdaQ) + crossprod(tcrossprod(A, LambdaQ))
  # tcrossprod(Q) + crossprod(tcrossprod(A, Q))
  # Matrix::rowSums(Q^2) + Matrix::colSums(tcrossprod(Q, A)^2)
}

get_V <- function(x) {
  # compute the variance of the outcome y (not including sigma^2)
  A <- lme4::getME(x, "A")
  I <- diag(nobs(x))
  crossprod(A) + I
}

#' Get the square root of a matrix using eigenvalue decomposition, and solve
#' for a linear system
#'
#' Solve for \eqn{x} in the linear system
#' \eqn{Ax = b}, where \eqn{A} is a
#' symmetric matrix square root of \eqn{M} with eigenvalue
#' decomposition such that \eqn{AA = M}.
#'
#' @param M a symmetric positive definite/semi-definite matrix
#' @param b a numeric vector or matrix
solve_eigen_sqrt <- function(M, b) {
  ei <- eigen(M)
  d <- ei$values
  dsqrtinv <- 1 / sqrt(d)
  dsqrtinv[d <= 0] <- 0
  V <- ei$vectors
  Msqrtinv <- V %*% diag(dsqrtinv) %*% t(V)
  Msqrtinv %*% b
}

#' @importFrom methods as
get_reflate_b <- function(x) {
  # Extract required quantities from the S4 object
  PR <- x@pp
  # rsp <- x@resp
  L <- PR$L()
  RX <- PR$RX()
  Lambdat <- PR$Lambdat
  A <- Lambdat %*% PR$Zt
  X <- model.matrix(x)
  u <- x@u

  Js <- lme4::ngrps(x)
  qs <- lengths(x@cnms)
  nqs <- Js * qs
  nqseq <- rep.int(seq_along(nqs), nqs)

  # Hat matrix for u
  Vb <- var_blup(L, RX, A, Lambdat, X)
  Vbstar <- Vb
  Vbstar[crossprod(Lambdat) == 0] <- 0
  Vbstar <- as(as(Vbstar, "symmetricMatrix"), "CsparseMatrix")
  R_Vbstar <- try(Matrix::chol(Vbstar), silent = TRUE)
  if (inherits(R_Vbstar, "try-error")) {
    bstar <- crossprod(Lambdat,
                       solve_eigen_sqrt(Vbstar, crossprod(Lambdat, u)))
  } else {
    bstar <- crossprod(Lambdat,
                       Matrix::solve(t(R_Vbstar), crossprod(Lambdat, u)))
  }
  bstar_lst <- split(bstar, nqseq)
  ml <- lapply(seq_along(bstar_lst),
               # easier to work with the transposed version
               function(i) matrix(bstar_lst[[i]], nrow = qs[i]))
  return(ml)
}

get_reflate_e <- function(x) {
  rsp <- x@resp
  hat_e <- hatvalues(x, fullHatMatrix = FALSE)
  (rsp$y - rsp$mu) / sqrt(1 - hat_e)
}

get_reflate_b_cgr <- function(x) {
  b <- lme4::getME(x, "b")
  Js <- lme4::ngrps(x)
  qs <- lengths(x@cnms)

  LR <- lme4::vec2mlist(x@theta, n = qs, symm = FALSE)

  # Hat matrix for u
  nqs <- Js * qs
  nqseq <- rep.int(seq_along(nqs), nqs)

  b_lst <- split(b, nqseq)
  ml <- lapply(seq_along(b_lst), function(i) {
    b_mat <- matrix(b_lst[[i]], nrow = qs[i])
    S <- tcrossprod(b_mat - rowMeans(b_mat)) / Js[i]
    LS <- try(t(Matrix::chol(S)), silent = TRUE)
    if (inherits(LS, "try-error")) {
      LR[[i]] %*% solve_eigen_sqrt(S, b_mat) * sigma(x)
    } else {
      LR[[i]] %*% Matrix::solve(LS, b_mat) * sigma(x)
    }
  })
  # }
  return(ml)
}

get_reflate_e_cgr <- function(x) {
  estar <- resid(x)
  estar / sd(estar) * sigma(x)
}

get_zeta <- function(r, R) {
  # R is the cholesky factor of V
  Zeta <- Matrix::solve(t(R), r)
  Zeta - mean(Zeta)
}

get_zeta_eigen <- function(r, V) {
  Zeta <- solve_eigen_sqrt(V, r)
  Zeta - mean(Zeta)
}
