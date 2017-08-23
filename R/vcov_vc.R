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

# theta_to_Lambdat <- function(theta, Js, qs) {
#   stopifnot(length(Js) == length(qs))
#   LR <- lme4::vec2mlist(x@theta, n = qs, symm = FALSE)
#   Ldt_lst <- lapply(seq_along(Js),
#                     function(i) Matrix::Diagonal(Js[[i]]) %x% LR[[i]])
#   t(Matrix::bdiag(Ldt_lst))
# }

#' Asymptotic Covariance Matrix for Random Effects
#'
#' Return the asymptotic covariance matrix of random effect standard
#' deviations (or variances) for a fitted model object, using the Hessian
#' evaluated at the (restricted) maximum likelihood estimates.
#'
#' Although it's easy to obtain the Hessian for \eqn{\theta}, the relative
#' Cholesky factor, in \pkg{lme4}, there is no easy way to obtain the Hessian
#' for the variance components. This function uses \code{\link{devfun_mer}()} to
#' obtain the Hessian (\eqn{H}) of variance components (or standard deviations,
#' SD), and then obtain the asymptotic covariance matrix as \eqn{-2 H^{-1}}.
#'
#' @param x A fitted merMod object from \code{\link[lme4]{lmer}}.
#' @param sd_cor Logical indicating whether to return asymptotic covariance
#'   matrix on SD scale (if \code{TRUE}) or on variance scale
#'   (if \code{FALSE}).
#' @param print_names Logical, whether to print the names for the covariance
#'   matrix.
#' @return A (q + 1) * (q + 1) symmetric matrix of the covariance
#'   matrix of (\eqn{\tau, \sigma}) (if \code{sd_cor = TRUE}) or
#'   (\eqn{\tau^2, \sigma^2}) (if \code{sd_cor = FALSE}), where q is the
#'   the number of estimated random-effects components (excluding \eqn{\sigma}).
#'   For example, for a model with random slope, \eqn{\tau} =
#'   (intercept SD, intercept-slope correlation, slope SD).
#' @seealso \code{\link[lme4]{vcov.merMod}} for covariance matrix of fixed
#'   effects, \code{\link[lme4]{confint.merMod}} for confidence intervals of all
#'   parameter estimates, and \code{\link{devfun_mer}} for the underlying
#'   function to produce the deviance function.
#' @export
#' @examples
#' library(lme4)
#' data(Orthodont, package = "nlme")
#' fm1 <- lmer(distance ~ age + (age | Subject), data = Orthodont)
#' vc <- VarCorr(fm1)
#' # Standard deviation only
#' print(vc, comp = c("Std.Dev"))
#' # Asymptotic variance-covariance matrix of (tau, sigma):
#' vcov_vc(fm1, sd_cor = TRUE)
#' # Compare with (parametric) bootstrap results :
#' get_sdcor <- function(x) {
#'   as.data.frame(lme4::VarCorr(x), order = "lower.tri")[ , "sdcor"]
#' }
#' boo <- bootstrap_mer(fm1, get_sdcor, type = "parametric", nsim = 200L)
#' # There might be failures in some resamples
#' cov(boo$t, use = "complete.obs")
vcov_vc <- function(x, sd_cor = TRUE, print_names = TRUE) {
  if (!lme4::isLMM(x)) {
    stop("currently only linear mixed model of class `merMod` is supported")
  }
  dd <- devfun_mer(x)
  n_th <- length(x@theta)
  qs <- lengths(x@cnms)
  if (sd_cor) {
    # from_chol <- lme4::Cv_to_Sv
    to_chol <- lme4::Sv_to_Cv
  } else {
    # from_chol <- lme4::Cv_to_Vv
    to_chol <- lme4::Vv_to_Cv
  }
  dd2 <- function(vc) {
    sigma <- if (sd_cor) vc[n_th + 1] else sqrt(vc[n_th + 1])
    th_sig <- c(to_chol(vc, n = qs, s = sigma), sigma)
    dd(th_sig)
  }
  vc <- lme4::VarCorr(x)
  vdd <- as.data.frame(vc, order = "lower.tri")
  vc_pars <- if (sd_cor) vdd[ , "sdcor"] else vdd[ , "vcov"]
  # vc <- from_chol(x@theta, n = qs, s = sigma(x))
  hess <- numDeriv::hessian(dd2, vc_pars)
  vv <- 2 * Matrix::solve(hess)
  if (print_names) {
    # prefix <- if (sd_cor) {
    #   c("sd_", "cor_", "sigma")
    # } else {
    #   c("var_", "cov_", "sigma2")
    # }
    # nms <- with(vdd,
    #             ifelse(!is.na(var1) & !is.na(var2),
    #                    paste(prefix[2], var2, ".", var1, "|", grp, sep = ""),
    #                    ifelse(grp == "Residual", prefix[3],
    #                           paste(prefix[1], var1, "|", grp, sep = ""))))
    nms <- make_vcnames(vdd, sd_cor = sd_cor)
    dimnames(vv) <- list(nms, nms)
  }
  vv
}

make_vcnames <- function(x, sd_cor = TRUE) {
  # If `x` is of class `VarCorr.merMod`
  if (inherits(x, "VarCorr.merMod")) {
    vdd <- as.data.frame(x)
  } else if (inherits(x, "data.frame")) {
    vdd <- x
  } else {
    stop("input object is not obtained from lme4::VarCorr()")
  }
  prefix <- if (sd_cor) {
    c("sd_", "cor_", "sigma")
  } else {
    c("var_", "cov_", "sigma2")
  }
  with(vdd,
       ifelse(!is.na(var1) & !is.na(var2),
              paste(prefix[2], var2, ".", var1, "|", grp, sep = ""),
              ifelse(grp == "Residual", prefix[3],
                     paste(prefix[1], var1, "|", grp, sep = ""))))
}

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

#' Asymptotic Covariance Matrix for Cholesky Factor of Random Effects
#'
#' @param x A fitted merMod object from \code{\link[lme4]{lmer}}.
#' @seealso \code{\link{vcov_vc}}
#' @importFrom numDeriv hessian
#' @importFrom stats update
#' @export
vcov_theta <- function(x) {
  org_data <- x@frame
  x_devfun <- update(x, data = org_data, devFunOnly = TRUE)
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
#' @examples
#' library(lme4)
#' fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
#' dd <- devfun_mer(fm01ML)
#' # Asymptotic variance-covariance matrix of (theta, sigma):
#' 2 * solve(numDeriv::hessian(dd, c(fm01ML@theta, sigma(fm01ML))))
devfun_mer <- function(x) {
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

#' @rdname devfun_mer
devfun_mer2 <- function(x) {
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

    function(theta) {
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
      logDet + df * (1 + log(2 * pi * pwrss) - log(df))
    }
  })
}

format_perc <- function (probs, digits = 3) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE,
               digits = digits), "%")
}

prof_ci_icc <- function(x, level = 0.95) {
  dd <- devfun_mer2(x)
  dd_na <- function(x) tryCatch(dd(x), error = function(e) NA)
  th0 <- x@theta
  min_dd <- dd(th0)
  fup <- function(eps) dd_na(th0 + eps) - min_dd - qchisq(level, 1)
  ul <- th0 + uniroot(fup, c(0, 1e6))$root
  flow <- function(eps) dd_na(th0 - eps) - min_dd - qchisq(level, 1)
  ll <- try(th0 - uniroot(flow, c(0, th0))$root, silent = TRUE)
  if (inherits(ll, "try-error")) {
    ll <- 0
  }
  out <- 1 / (1 + c(ll, ul)^(-2))
  a <- .5 + c(-1, 1) * level / 2
  names(out) <- format_perc(a)
  out
}
