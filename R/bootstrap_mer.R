#' Run Various Bootstrap for Mixed Models.
#'
#' Run multilevel parametric, residual, and case bootstrap with different
#'   options
#'
#' \code{bootstrap_mer} performs different bootstrapping methods to fitted
#' model objects using the lme4 package. Currently, only models fitted
#' using \code{\link[lme4]{lmer}} is supported.
#' @param x A fitted \code{merMod} object from \code{\link[lme4]{lmer}}.
#' @param FUN A function taking a fitted \code{merMod} object as input and
#'   returning the statistic of interest, which must be a (possibly named)
#'   numeric vector.
#' @param nsim A positive integer telling the number of simulations, positive
#'   integer; the bootstrap \eqn{R}.
#' @param seed Optional argument to \code{\link[base]{set.seed}}.
#' @param type A character string indicating the type of multilevel bootstrap.
#'   Currently, possible values are \code{"parametric"}, \code{"residual"},
#'   \code{"residual_cgr"}, \code{"residual_trans"}, or \code{"case"}.
#' @param lv1_resample Logical indicating whether to sample with replacement
#'   the level-1 units for each level-2 cluster. Only used for
#'   \code{type = "case"}. Default is \code{FALSE}.
#' @param corrected_trans Logical indicating whether to use the correct
#'   variance-covariance matrix of the residuals. If \code{FALSE}, use the
#'   variance of \eqn{y}; if \code{TRUE}, use the variance of \eqn{y - X \hat
#'   \beta}. Only used for \code{type = "residual_trans"}.
#' @param .progress Logical indicating whether to display progress bar (using
#'   \code{\link[utils]{txtProgressBar}}).
#' @param verbose Logical indicating if progress should print output.
#' @return An object of S3 class \code{"boot"}, compatible with \pkg{boot}
#'   package's \code{\link[boot]{boot}()}. It contains the following components:
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
#' @importFrom lme4 refit lmer lmerControl
#' @importFrom stats formula
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' library(lme4)
#' fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
#' mySumm <- function(x) {
#'   c(getME(x, "beta"), sigma(x))
#' }
#' # Covariance preserving residual bootstrap
#' boo01 <- bootstrap_mer(fm01ML, mySumm, type = "residual", nsim = 100)
#' # Plot bootstrap distribution of fixed effect
#' library(boot)
#' plot(boo01, index = 1)
#' # Get confidence interval
#' boot.ci(boo01, index = 2, type = c("norm", "basic", "perc"))
#' # BCa using influence values computed from `empinf_`
#' boot.ci(boo01, index = 2, type = "bca", L = empinf_mer(fm01ML, mySumm, 2))
bootstrap_mer <- function(x, FUN, nsim = 1, seed = NULL,
                          type = c("parametric", "residual", "residual_cgr",
                                   "residual_trans", "case"),
                          corrected_trans = FALSE,
                          lv1_resample = FALSE, .progress = FALSE,
                          verbose = FALSE) {
  type <- match.arg(type)
  if (type == "parametric") {
    return(lme4::bootMer(x, FUN, nsim, seed = seed, use.u = FALSE,
                         type = "parametric", verbose = FALSE))
  } else {
    if (!lme4::isLMM(x)) {
      stop("currently only linear mixed model of class `merMod` is supported")
    }
    # if (!identical(x@optinfo$conv$lme4, list())) {
    #   stop("The original model has convergence issue")
    # }
    stopifnot( (nsim <- as.integer(nsim[1])) > 0)
    if (.progress) {
      pb <- txtProgressBar(style = 3)
    }
    FUN <- match.fun(FUN)
    # if (!is.null(seed)) set.seed(seed)
    if (!missing(seed)) set.seed(seed)

    t0 <- FUN(x)
    if (!is.numeric(t0)) {
      stop("currently only handles functions that return numeric vectors")
    }

    # can use the switch function
    if (type %in% c("residual", "residual_cgr", "residual_trans")) {
      mle <- list(beta = x@beta, theta = x@theta)
      out <- list(sim = "parametric", ran.gen = NULL, mle = mle)
      ss <- .resid_resample(x, nsim, type = type, corrected = corrected_trans)
      ffun <- local({
        FUN
        refit
        x
        ss
        verbose
        length_t0 <- length(t0)
        function(i) {
          ret <- tryCatch(FUN(refit(x, ss[[i]])),
                          error = function(e) e)
          if (verbose) {
            cat(sprintf("%5d :", i))
            utils::str(ret)
          }
          if (.progress) {
            setTxtProgressBar(pb, i / nsim)
          }
          if (inherits(ret, "error"))
            structure(rep(NA, length_t0), fail.msgs = ret$message)
          else ret
        }
      })
    } else if (type == "case") {
      out <- list(sim = "ordinary", strata = rep(1, nobs(x)))
      ss <- .case_resample(x, nsim, lv1_resample = lv1_resample)
      ffun <- local({
        FUN
        formula_x <- formula(x)
        ss
        verbose
        length_t0 <- length(t0)
        use_REML <- as.logical(lme4::getME(x, "REML"))
        function(i) {
          df_i <- ss[[i]]
          ret <- tryCatch(
            FUN(lmer(formula_x, data = df_i, REML = use_REML,
                     control = lmerControl(calc.derivs = FALSE))),
            error = function(e) e)
          if (verbose) {
            cat(sprintf("%5d :", i))
            utils::str(ret)
          }
          if (.progress) {
            setTxtProgressBar(pb, i / nsim)
          }
          if (inherits(ret, "error"))
            structure(rep(NA, length_t0), fail.msgs = ret$message)
          else ret
        }
      })
    }

    res <- lapply(seq_along(ss), ffun)
    # t.star <- matrix(unlist(res), nsim, length(t0), byrow = TRUE)
    # colnames(t.star) <- names(t0)
    t.star <- do.call(rbind, res)

    # Number of failed bootstrap

    boo <- structure(c(list(t0 = t0, t = t.star, R = nsim, data = x@frame,
                            seed = .Random.seed, statistic = FUN,
                            call = match.call()), out),
                     class = "boot")
  }
  return(boo)
}

# pupcross <- haven::read_sas(
#   "https://stats.idre.ucla.edu/wp-content/uploads/2016/02/pupcross.sas7bdat")
#
# library(lme4)
# x <- lmer(ACHIEV ~ PUPSEX + PUPSES + (PUPSES | PSCHOOL) + (1 | SSCHOOL),
#           data = pupcross)
# boo1 <- bootstrap_mer(x, function(x) x@theta * sigma(x), type = "resid",
#                       nsim = 1000)
# boo2 <- bootstrap_mer(x, function(x) x@theta * sigma(x), type = "resid_trans",
#                       nsim = 1000)
# boo3 <- bootstrap_mer(x, function(x) x@theta * sigma(x), type = "resid_cgr",
#                       nsim = 1000)
