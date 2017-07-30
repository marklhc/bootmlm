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
  if (all(u == 0)) {  # may break down when variance is 0 for just one component
    bstar <- u
  } else {
    Vb <- var_blup(L, RX, A, Lambdat, X)
    Vbstar <- Vb
    Vbstar[crossprod(Lambdat) == 0] <- 0
    Vbstar <- as(as(Vbstar, "symmetricMatrix"), "CsparseMatrix")
    R_Vbstar <- Matrix::chol(Vbstar)
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
  if (all(b == 0)) {  # may break down when variance is 0 for just one component
    ml <- lapply(seq_along(qs), function(i) {
      matrix(0, nrow = qs[i], ncol = Js[i])
    })
  } else {
    nqs <- Js * qs
    nqseq <- rep.int(seq_along(nqs), nqs)

    b_lst <- split(b, nqseq)
    ml <- lapply(seq_along(b_lst), function(i) {
      b_mat <- matrix(b_lst[[i]], nrow = qs[i])
      LS <- t(Matrix::chol(tcrossprod(b_mat - rowMeans(b_mat)) / Js[i]))
      LR[[i]] %*% Matrix::solve(LS, b_mat) * sigma(x)
    })
  }
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
