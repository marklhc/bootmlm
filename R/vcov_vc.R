# vcov_vc <- function(x) {
#   # http://rstudio-pubs-static.s3.amazonaws.com/28864_dd1f084207d54f5ea67c9d1a9c845d01.html
#   if (isREML(x)) {
#     warning("refitting model with ML")
#     x <- refitML(x)
#   }
#   # if (!require("numDeriv")) stop("numDeriv package required")
#   vc <- lme4::VarCorr(x)
#   useSc <- attr(vc, "useSc")
#   dd <- lme4:::devfun2(x, useSc = useSc, signames = FALSE)
#   vdd <- as.data.frame(vc, order = "lower.tri")
#   pars <- vdd[, "sdcor"]
#   npar0 <- length(pars)
#   if (isGLMM(x)) {
#     pars <- c(pars, fixef(x))
#     hh1 <- numDeriv::hessian(dd, pars)
#     vv2 <- 2 * solve(hh1)
#     vv2 <- vv2[1:npar0, 1:npar0, drop = FALSE]
#   } else {
#     hh1 <- numDeriv::hessian(dd, pars)
#     vv2 <- 2 * solve(hh1)
#   }
#   nms <- apply(vdd[ , 1:3], 1,
#                function(x) paste(na.omit(x), collapse = "."))
#   dimnames(vv2) <- list(nms, nms)
#   return(vv2)
# }

devfun_sig <- function(sigma, .x) {
  sigsq <- sigma^2
  if (lme4::isREML(.x)) {
    df <- nobs(.x)
  } else {
    df <- nobs(.x) - length(.x@beta)
  }
  (.x@resp$wrss() + .x@pp$sqrL(1)) / sigsq + df *
    log(2 * pi * sigsq)
}

#' @importFrom numDeriv hessian
#' @importFrom stats update
vcov_theta <- function(x) {
  x_devfun <- update(x, devFunOnly = TRUE)
  hess <- numDeriv::hessian(x_devfun, x@theta)
  2 * Matrix::solve(hess)
}

# theta_to_Lambdat <- function(theta, Js, qs) {
#   stopifnot(length(Js) == length(qs))
#   LR <- lme4::vec2mlist(x@theta, n = qs, symm = FALSE)
#   Ldt_lst <- lapply(seq_along(Js),
#                     function(i) Matrix::Diagonal(Js[[i]]) %x% LR[[i]])
#   t(Matrix::bdiag(Ldt_lst))
# }

#' Deviance Function for Multilevel Models
#'
#' Creates a deviance function for a fitted model object, using \eqn{\theta}
#' and \eqn{sigma} as the parameters.
#'
#' The built-in function(s) in \pkg{lme4} for generating deviance function are
#' not exported and rely on C++ code. This function is mainly used to obtain
#' Hessian and asymptotic covariance matrix of the random effects.
#'
#' @param x A fitted merMod object from \code{\link[lme4]{lmer}}.
#' @return A deviance function with one argument \code{th_sig} that takes input
#'   of a numeric vector corresponding to the some estimated values of
#'   \eqn{\theta} and \eqn{\sigma}
#' @references Bates, D., M\"{a}chler, M., Bolker, B. M., & Walker, S. C.
#'   Fitting linear mixed-effects models using lme4. Retrieved from
#'   \url{https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf}
#' @references Implementation in \pkg{lme4pureR}:
#'   \url{https://github.com/lme4/lme4pureR/blob/master/R/pls.R}
#' @export
devfun2_mer <- function(x) {
  res <- x@resp
  offset <- res$offset
  y <- res$y
  PR <- x@pp
  X <- PR$X
  Zt <- PR$Zt
  n <- length(y)
  n_th <- length(x@theta)

  sqrtW <- Matrix::Diagonal(x = res$weights)
  WX <- sqrtW %*% X
  Wy <- sqrtW %*% y
  ZtW <- Zt %*% sqrtW
  XtWX <- crossprod(WX)
  XtWy <- crossprod(WX, Wy)
  ZtWX <- ZtW %*% WX
  ZtWy <- ZtW %*% Wy

  REML <- lme4::isREML(x)
  Lind <- PR$Lind

  local({  # mutable values stored in local environment
    Lambdat <- PR$Lambdat
    L <- Matrix::Cholesky(tcrossprod(Lambdat %*% ZtW), LDL = FALSE, Imult = 1)
    df <- n

    function(th_sig) {
      theta <- th_sig[1:n_th]
      sigma <- th_sig[n_th + 1]
      Lambdat@x <- theta[Lind]
      L <- Matrix::update(L, Lambdat %*% ZtW, mult = 1)
      cu <- Matrix::solve(L, Matrix::solve(L, Lambdat %*% ZtWy, system = "P"),
                          system = "L")
      RZX <- Matrix::solve(L, Matrix::solve(L, Lambdat %*% ZtWX, system = "P"),
                           system = "L")
      RXtRX <- as(XtWX - crossprod(RZX), "dpoMatrix")
      betahat <- Matrix::solve(RXtRX, XtWy - crossprod(RZX, cu))
      u <- Matrix::solve(L,
                         Matrix::solve(L, cu - RZX %*% betahat, system = "Lt"),
                         system = "Pt")
      b <- crossprod(Lambdat, u)
      mu <- crossprod(Zt, b) + X %*% betahat + offset
      wtres <- sqrtW %*% (y - mu)
      pwrss <- sum(wtres^2) + sum(u^2)
      logDet <- 2 * Matrix::determinant(L, logarithm = TRUE)$modulus
      if (REML) {
        logDet <- logDet + Matrix::determinant(RXtRX, logarithm = TRUE)$modulus
        df <- df - length(betahat)
      }
      attributes(logDet) <- NULL
      sigsq <- sigma^2
      logDet + df * log(2 * pi * sigsq) + pwrss / sigsq
    }
  })
}

# Example:
# dd2 <- devfun2_mer(m2)
# 2 * solve(numDeriv::hessian(dd2, c(m2@theta, sigma(m2))))
