#' Bootstrap confidence intervals for Two-Level Mixed Models
#'
#' This is a wrapper for getting CIs for multiple parameters after running
#' \code{\link{bootstrap_mer}}, to be consistent with a similar method for the
#' \code{\link[lme4]{bootMber}} class.
#'
#' @param object an object returned by \code{\link{bootstrap_mer}}.
#' @param parm a specification of which parameters are to be given confidence
#'             intervals as a vector of numbers. If missing, all parameters
#'             are considered.
#' @param level the confidence level required.
#' @param type character indicating the type of intervals required, as described
#'             in \code{\link[boot]{boot.ci}}. Currently \code{"stud"} is not
#'             supported.
#' @param L empirical influence values required for \code{type = "bca"} as
#'          described in \code{\link[boot]{boot.ci}}.
#' @param ... additional argument(s) passed to \code{\link[boot]{boot.ci}}.
#' @importFrom stats confint
#' @return A matrix (or vector) with columns giving lower and upper confidence
#'         limits for each parameter.
#' @export
#' @examples
#' library(lme4)
#' fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
#' mySumm <- function(x) {
#'   c(getME(x, "beta"), sigma(x))
#' }
#'
#' # residual bootstrap
#' boo_resid <- bootstrap_mer(fm01ML, mySumm, type = "residual", nsim = 100)
#' confint(boo_resid, type = "bca", L = empinf_merm(fm01ML, mySumm))

# parametric bootstrap
boo_parametric <- bootstrap_mer(fm01ML, mySumm, type = "parametric", nsim = 100)
class(boo_parametric) # returns "bootMer" "boot"
confint(boo_parametric, type = "norm") # returns a confidence interval

# residual bootstrap
boo_resid <- bootstrap_mer(fm01ML, mySumm, type = "residual", nsim = 100)
class(boo_resid) # returns "boot"
confint(boo_resid, type = "norm") # does not work
boot.ci(boo_resid, type = "norm") # kind of works but only for one value at a time
confint.boot <- function(object, parm,
                         level = 0.95,
                         type = c("norm", "basic", "perc", "bca"),
                         L = NULL,
                         ...) {
  # Need to consider how to get `stud` to work
  type <- match.arg(type)
  if (missing(parm)) {
    parm <- seq_along(object$t0)
  }
  bnm <- switch(type, norm = "normal", basic = "basic", perc = "percent",
                bca = "bca")
  bind <- if (type == "norm") 2:3 else 4:5
  if (type == "bca" && is.null(L)) {
    stop("BCa CIs require inputs of influence values. ",
         "The `empinf_mer()` function can be used to estimate those using the ",
         "grouped jackknife.")
  }
  btab <- t(vapply(parm, function(i) {
    boot::boot.ci(object, index = i, conf = level, type = type,
                  L = L, ...)[[bnm]][bind]
    }, FUN.VALUE = numeric(2)))
  rownames(btab) <- names(object$t0)
  colnames(btab) <-
    paste(format(100 * c(1 - level, 1 + level) / 2,
                 trim = TRUE, scientific = FALSE,
                 digits = 3), "%")
  return(btab)
}
