format_perc <- function (probs, digits = 3) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE,
               digits = digits), "%")
}

#' Profile Likelihood Confidence Interval for Intraclass Correlation
#'
#' Compute confidence intervals for the intraclass correlation of a model
#' fit of class \code{\link[lme4]{merMod-class}}.
#'
#' This function uses the \code{\link[stats]{uniroot}} function and determine
#' the lower and upper limit by evaluating the profile deviance (for ML) or the
#' -2 profile REML criterion. It works by obtaining the interval for
#' \eqn{\theta = \tau / \sigma} and transforming the two limits.
#'
#' The resulting interval is bounded by zero for its lower limit.
#'
#' @param x A fitted merMod object from \code{\link[lme4]{lmer}}.
#' @param level Confidence level between 0 and 1. Default is .95.
#' @param eps_max The maximum value that the upper limit can be. Default is
#'   1e6.
#' @return A named numeric vector with two values showing the lower and upper
#'   limit of the intraclass correlation.
prof_ci_icc <- function(x, level = 0.95, eps_max = 1e6) {
  dd <- devfun_mer2(x)
  dd_na <- function(x) tryCatch(dd(x), error = function(e) NA)
  th0 <- x@theta
  min_dd <- dd(th0)
  fup <- function(eps) dd_na(th0 + eps) - min_dd - stats::qchisq(level, 1)
  ul <- th0 + stats::uniroot(fup, c(0, eps_max))$root
  flow <- function(eps) dd_na(th0 - eps) - min_dd - stats::qchisq(level, 1)
  ll <- try(th0 - stats::uniroot(flow, c(0, th0))$root, silent = TRUE)
  if (inherits(ll, "try-error")) {
    ll <- 0
  }
  out <- 1 / (1 + c(ll, ul)^(-2))
  a <- .5 + c(-1, 1) * level / 2
  names(out) <- format_perc(a)
  out
}
