#' Run Various Bootstrap for Mixed Models.
#'
#' Run multilevel parametric, residual, and case bootstrap with different
#'   options
#'
#' \code{bootstrap_mer} performs different bootstrapping methods to fitted
#' model objects using the \pkg{lme4} package. Currently, only models fitted
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
#'   \code{"residual_cgr"}, \code{"residual_trans"}, \code{"reb"},
#'   or \code{"case"}.
#' @param lv1_resample Logical indicating whether to sample with replacement
#'   the level-1 units for each level-2 cluster. Only used for
#'   \code{type = "case"}. Default is \code{FALSE}.
#' @param corrected_trans Logical indicating whether to use the correct
#'   variance-covariance matrix of the residuals. If \code{FALSE}, use the
#'   variance of \eqn{y}; if \code{TRUE}, use the variance of \eqn{y - X \hat
#'   \beta}. Only used for \code{type = "residual_trans"}.
#' @param reb_scale Logical indicating whether to scale the residuals for the
#'   random effect block bootstrap
#' @param .progress Logical indicating whether to display progress bar (using
#'   \code{\link[utils]{txtProgressBar}}).
#' @param verbose Logical indicating if progress should print output.
#' @return An object of S3 class \code{"boot"}, compatible with \pkg{boot}
#'   package's \code{\link[boot]{boot}()}. It contains the following components:
#'
#'   \item{t0}{The original statistic from \code{FUN(x)}.}
#'   \item{t}{A matrix with \code{nsim} rows containing the bootstrap
#'   distribution of the statistic.}
#'   \item{R}{The value of \code{nsim} passed to the function.}
#'   \item{data}{The data used in the original analysis.}
#'   \item{seed}{The value of \code{.Random.seed} when \code{bootstrap_mer}
#'   started to work.}
#'   \item{statistic}{The function \code{FUN} passed to \code{bootstrap_mer}.}
#'
#' See the documentation in for \code{link[boot]{boot}()} for the other
#'   components.
#' @references Carpenter, J. R., Goldstein, H., & Rasbash, J. (2003). A novel
#'   bootstrap procedure for assessing the relationship between class size and
#'   achievement. Journal of the Royal Statistical Society. Series C (Applied
#'   Statistics), 52, 431–443. https://doi.org/10.1111/1467-9876.00415
#' @references Chambers, R., & Chandra, H. (2013). A random effect block
#'   bootstrap for clustered data. Journal of Computational and Graphical
#'   Statistics, 22(2), 452–470. https://doi.org/10.1080/10618600.2012.681216
#' @references Davison, A. C. and Hinkley, D. V. (1997). Bootstrap methods and
#'   their application. Cambridge, UK: Cambridge University Press.
#' @references Morris, J. S. (2002). The BLUPs are not "best" when it comes to
#'   bootstrapping. Statistics & Probability Letters, 56(4), 425–430.
#'   https://doi.org/10.1016/S0167-7152(02)00041-X
#' @references Van der Leeden, R., Meijer, E., & Busing, F. M. T. A. (2008).
#'   Resampling multilevel models. In J. de Leeuw & E. Meijer (Eds.), Handbook
#'   of multilevel Analysis (pp. 401–433). New York, NY: Springer.
#' @seealso
#' \itemize{
#'   \item \code{\link[boot]{boot}} for single-level bootstrapping,
#'   \item \code{\link[lme4]{bootMer}} for parametric and semi-parametric
#'   bootstrap implemented in lme4, and
#'   \item \code{\link[boot]{boot.ci}} for getting bootstrap confidence
#'   intervals and \code{\link[boot]{plot.boot}} for plotting the bootstrap
#'   distribution.}
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
#' # BCa using influence values computed from `empinf_mer`
#' boot.ci(boo01, index = 2, type = "bca", L = empinf_mer(fm01ML, mySumm, 2))
bootstrap_mer <- function(x, FUN, nsim = 1, seed = NULL,
                          type = c("parametric", "residual", "residual_cgr",
                                   "residual_trans", "reb", "case"),
                          corrected_trans = FALSE,
                          lv1_resample = FALSE,
                          reb_scale = FALSE,
                          .progress = FALSE, verbose = FALSE) {
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
    if (type %in% c("residual", "residual_cgr", "residual_trans", "reb")) {
      mle <- list(beta = x@beta, theta = x@theta)
      out <- list(sim = "ordinary", ran.gen = NULL, mle = mle)
      if (type == "reb") {
        ss <- .reb_resample(x, nsim, scale = reb_scale)
      } else {
        ss <- .resid_resample(x, nsim, type = type, corrected = corrected_trans)
      }
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
        # x
        FUN
        formula_x <- formula(x)
        # update
        ss
        verbose
        length_t0 <- length(t0)
        use_REML <- lme4::isREML(x)
        function(i) {
          df_i <- ss[[i]]
          ret <- tryCatch(
            FUN(lmer(formula_x, data = df_i, REML = use_REML,
                     control = lmerControl(calc.derivs = FALSE))),
            # FUN(update(x, data = df_i,
            #            control = lmerControl(calc.derivs = FALSE))),
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
  attr(boo, "boot_type") <- "boot"
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
