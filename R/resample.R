#' @importFrom Matrix crossprod tcrossprod t mean
#' @importFrom stats nobs resid runif sd sigma hatvalues model.matrix
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
  Zt <- PR$Zt
  X <- PR$X
  fixed <- X %*% x@beta

  Js <- lme4::ngrps(x)
  qs <- lengths(x@cnms)
  nqs <- Js * qs
  nqseq <- rep.int(seq_along(nqs), nqs)

  bstar <- get_reflate_b(x)
  bstar_lst <- split(bstar, nqseq)
  ml <- lapply(seq_along(bstar_lst),
               # easier to work with the transposed version
               function(i) matrix(bstar_lst[[i]], nrow = qs[i]))

  estar <- get_reflate_e(x)

  replicate(nsim, {
    # faster!, and easier to read
    bstar_new <- unlist(lapply(ml, function(m) {
      m[ , sample.int(ncol(m), replace = TRUE)]
    }))

    as.vector(fixed + crossprod(Zt, bstar_new) + sample(estar, replace = TRUE))
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
  X <- unname(PR$X)
  fixed <- X %*% x@beta

  ml <- get_reflate_b_cgr(x)
  estar <- get_reflate_e_cgr(x)

  replicate(nsim, {
    bstar_new <- unlist(lapply(ml, function(m) {
      m[ , sample.int(ncol(m), replace = TRUE)]
    }))

    as.vector(fixed + crossprod(Zt, bstar_new) + sample(estar, replace = TRUE))
  },
  simplify = FALSE)
}

.resid_trans_resample <- function(x, nsim = 1, seed = NULL,
                                  corrected = FALSE) {
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
  X <- PR$X
  fixed <- X %*% x@beta

  # Transform residuals
  r <- x@resp$y - fixed  # the variance of r may not be V
  V <- get_V(x)
  # var_r <- (V - X %*% Matrix::chol2inv(RX) %*% t(X)) * (V != 0)
  R <- Matrix::chol(V)
  if (corrected) {
    RX <- PR$RX()
    Vr <- (V - X %*% Matrix::chol2inv(RX) %*% t(X)) * (V != 0)
    Rr <- Matrix::chol(Vr)
    Zeta <- get_zeta(r, Vr)
  } else {
    Zeta <- get_zeta(r, R)  # can use a corrected version
  }

  replicate(nsim, {
    as.vector(fixed + crossprod(R, sample(Zeta, replace = TRUE)))
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
