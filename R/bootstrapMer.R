#' Run Various Bootstrap for Mixed Models.
#'
#' Run multilevel bootstrap with three options
#'
#' \code{bootstrapMer} performs different bootstrapping methods to fitted
#' model objects using the lme4 package. Currently, only models with fitted
#' using \code{\link[lme4]{lmer}} is supported.
#' @param x A fitted merMod object from lmer.
#' @param FUN A function taking a fitted merMod object as input and returning
#'   the statistic of interest, which must be a (possibly named) numeric vector.
#' @param nsim Number of simulations, positive integer; the bootstrap B (or R).
#' @param seed Optional argument to set.seed.
#' @param type A character string indicating the type of multilevel bootstrap.
#'   Possible values are "resid", "resid_cgr", or "resid_trans".
#' @param verbose Logical indicating if progress should print output
#' @return An object of S3 class "boot", compatible with boot package's boot()
#'   result. It contains the following components:
#'
#'   \item{t0}{The original statistic from \code{FUN(x)}.}
#'   \item{t}{A matrix with \code{nsim} rows containing the bootstrap
#'            distribution of the statistic.}
#'   \item{R}{The value of \code{nsim} passed to the function.}
#'   \item{data}{The data used in the original analysis.}
#' @seealso \code{\link[boot]{boot}} for single-level bootstrapping,
#'   \code{\link[lme4]{bootMer}} for parametric and semi-parametric bootstrap
#'   implemented in lme4, and \code{\link[boot]{boot.ci}} for getting
#'   bootstrap confidence intervals.
#' @importFrom lme4 refit
#' @export
#' @examples
#' library(lme4)
#' fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
#' # mathematically correct residual bootstrap
#' mySumm <- function(x) {
#'   c(getME(x, "beta"), sigma(x))
#' }
#' bootstrapMer(fm01ML, mySumm, type = "resid_trans", nsim = 100)
bootstrapMer <- function(x, FUN, nsim = 1, seed = NULL,
                         type = c("resid", "resid_cgr", "resid_trans"),
                         verbose = FALSE) {
  if (!lme4::isLMM(x)) {
    stop("currently only linear mixed model of class `merMod` is supported")
  }
  # if (!identical(x@optinfo$conv$lme4, list())) {
  #   stop("The original model has convergence issue")
  # }
  stopifnot( (nsim <- as.integer(nsim[1])) > 0)
  FUN <- match.fun(FUN)
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  # method <- match.arg(method)

  t0 <- FUN(x)
  if (!is.numeric(t0)) {
    stop("currently only handles functions that return numeric vectors")
  }

  mle <- list(beta = x@beta, theta = x@theta)

  # can use the switch function
  if (type == "resid") {
    ss <- .resid_resample(x, nsim)
  } else if (type == "resid_trans") {
    ss <- .resid_trans_resample(x, nsim)
  } else if (type == "resid_cgr") {
    ss <- .resid_cgr_resample(x, nsim)
  }

  ffun <- local({
    FUN
    refit
    x
    ss
    verbose
    length.t0 <- length(t0)
    function(i) {
      # ret <- tryCatch(FUN(refit.merMod2(x, ss[[i]],
      #                                   control = lmerControl(
      #                                     calc.derivs = calc.derivs))),
      #                 error = function(e) e)
      ret <- tryCatch(FUN(refit(x, ss[[i]])),
                      error = function(e) e)
      if (verbose) {
        cat(sprintf("%5d :", i))
        utils::str(ret)
      }
      if (inherits(ret, "error"))
        structure(rep(NA, length.t0), fail.msgs = ret$message)
      else ret
    }
  })
  res <- lapply(seq_along(ss), ffun)
  # t.star <- matrix(unlist(res), nsim, length(t0), byrow = TRUE)
  # colnames(t.star) <- names(t0)
  t.star <- do.call(rbind, res)

  # Number of failed bootstrap

  structure(list(t0 = t0, t = t.star, R = nsim, data = x@frame,
                 seed = .Random.seed, statistic = FUN,
                 sim = "parameteric", call = match.call(), ran.gen = NULL,
                 mle = mle), class = "boot")
}

# pupcross <- haven::read_sas(
#   "https://stats.idre.ucla.edu/wp-content/uploads/2016/02/pupcross.sas7bdat")
#
# library(lme4)
# x <- lmer(ACHIEV ~ PUPSEX + PUPSES + (PUPSES | PSCHOOL) + (1 | SSCHOOL),
#           data = pupcross)
# boo1 <- bootstrapMer(x, function(x) x@theta * sigma(x), type = "resid",
#                      nsim = 1000)
# boo2 <- bootstrapMer(x, function(x) x@theta * sigma(x), type = "resid_trans",
#                      nsim = 1000)
# boo3 <- bootstrapMer(x, function(x) x@theta * sigma(x), type = "resid_cgr",
#                      nsim = 1000)
