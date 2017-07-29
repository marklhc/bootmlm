#' @importFrom Matrix crossprod tcrossprod t mean
#' @importFrom stats nobs resid runif sd sigma
.resid_resample <- function(x, nsim = 1, seed = NULL) {
  # residual bootstrap
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1)
  }
  # RNGstate <- .Random.seed

  # Extract required quantities from the S4 object
  PR <- x@pp
  rsp <- x@resp
  L <- PR$L()
  RX <- PR$RX()
  RZX <- PR$RZX
  Lambdat <- PR$Lambdat
  # Zt <- PR$Zt
  # A <- unname(Lambdat %*% Zt)
  A <- unname(Lambdat %*% PR$Zt)
  X <- unname(PR$X)
  fixed <- X %*% x@beta
  u <- x@u

  Js <- lme4::ngrps(x)
  # Jseq <- seq_along(Js)
  qs <- lengths(x@cnms)
  nqs <- Js * qs
  nqseq <- rep.int(seq_along(nqs), nqs)

  # Hat matrix for u
  if (all(u == 0)) {  # may break down when variance is 0 for just one component
    ustar <- u
  } else {
    I <- diag(sum(nqs))
    Vu <- I - tcrossprod(
      Matrix::solve(L, cbind(I, t(Matrix::solve(t(RX), t(RZX)))),
                    system = "Lt")
    )
    ustar <- u / sqrt(Matrix::diag(Vu))
  }
  ustar_lst <- split(ustar, nqseq)
  ml <- lapply(seq_along(ustar_lst),
               # easier to work with the transposed version
               function(i) matrix(ustar_lst[[i]], nrow = qs[i]))

  # Hat matrix for e
  CL <- Matrix::solve(L, A, system = "L")
  CR <- Matrix::solve(t(RX), t(X) - crossprod(RZX, CL))
  hatvalues <- Matrix::colSums(CR^2) + Matrix::colSums(CL^2)
  estar <- (rsp$y - rsp$mu) / sqrt(1 - hatvalues)

  replicate(nsim, {
    # faster!, and easier to read
    ustar_new <- unlist(lapply(ml, function(m) {
      m[ , sample.int(ncol(m), replace = TRUE)]
    }))

    as.vector(
      fixed + crossprod(A, ustar_new) + sample(estar, replace = TRUE))
  },
  simplify = FALSE)
}

.resid_cgr_resample <- function(x, nsim = 1, seed = NULL) {
  # reflate residuals
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1)
  }
  # RNGstate <- .Random.seed

  # Extract required quantities from the S4 object
  PR <- x@pp

  Zt <- PR$Zt
  Lambdat <- PR$Lambdat
  X <- unname(PR$X)
  fixed <- X %*% x@beta
  u <- x@u
  b <- Matrix::crossprod(Lambdat, u)

  Js <- lme4::ngrps(x)
  qs <- lengths(x@cnms)
  sigma_x <- sigma(x)

  LR <- lme4::vec2mlist(x@theta, n = qs, symm = FALSE)

  # Hat matrix for u
  if (all(u == 0)) {  # may break down when variance is 0 for just one component
    ml <- lapply(seq_along(b_lst), function(i) {
      matrix(0, nrow = qs[i], ncol = Js[i])
    })
  } else {
    nqs <- Js * qs
    nqseq <- rep.int(seq_along(nqs), nqs)

    b_lst <- split(b, nqseq)
    ml <- lapply(seq_along(b_lst), function(i) {
      b_mat <- matrix(b_lst[[i]], nrow = qs[i])
      LS <- t(Matrix::chol(tcrossprod(b_mat - rowMeans(b_mat)) / Js[i]))
      LR[[i]] %*% Matrix::solve(LS, b_mat) * sigma_x
    })
  }
  estar <- resid(x)
  estar <- estar / sd(estar) * sigma_x

  replicate(nsim, {
    bstar_new <- unlist(lapply(ml, function(m) {
      m[ , sample.int(ncol(m), replace = TRUE)]
    }))

    as.vector(fixed + crossprod(Zt, bstar_new) + sample(estar, replace = TRUE))
  },
  simplify = FALSE)
}

.resid_trans_resample <- function(x, nsim = 1, seed = NULL) {
  # transform residuals to be independent
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1)
  }
  # RNGstate <- .Random.seed

  # Extract required quantities from the S4 object
  PR <- x@pp
  RX <- PR$RX()
  A <- unname(PR$Lambdat %*% PR$Zt)
  X <- unname(PR$X)
  fixed <- X %*% x@beta
  nobs_x <- nobs(x)

  # Transform residuals
  r <- x@resp$y - fixed  # the variance of r may not be V
  I <- diag(nobs_x)
  V <- crossprod(A) + I
  var_r <- (V - X %*% Matrix::chol2inv(RX) %*% t(X)) * (V != 0)
  R <- Matrix::chol(var_r)
  Zeta <- Matrix::solve(t(R), r)
  Zeta_c <- Zeta - mean(Zeta)

  replicate(nsim, {
    as.vector(fixed + crossprod(R, sample(Zeta_c, replace = TRUE)))
  },
  simplify = FALSE)
}

# pupcross <- haven::read_sas(
#   "https://stats.idre.ucla.edu/wp-content/uploads/2016/02/pupcross.sas7bdat")
#
# library(lme4)
# x <- lmer(ACHIEV ~ PUPSEX + PUPSES + (PUPSES | PSCHOOL) + (1 | SSCHOOL),
#           data = pupcross)
# .resid_resample(x, 1000)
# .resid_cgr_resample(x, 1000)
# .resid_trans_resample(x, 1000)
