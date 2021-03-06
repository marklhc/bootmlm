#' Empirical Influence Values for Two-Level Mixed Models
#'
#' This function calculates the empirical influence values for a statistic in a
#' given fitted model object using the delete-\eqn{m_j} jackknife.
#'
#' \code{empinf_mer} computes non-parametric influence function of models
#' fitted using \code{\link[lme4]{lmer}} by deleting one cluster at a time. See
#' van der Leeden, Meijer, and Busing (2008, pp. 420--422) for more information.
#' Whereas \code{empinf_mer} computes influence values for a specified position
#' (as specified with the \code{index} argument) of the output of \code{FUN},
#' \code{empinf_merm} computes influence values for every element in
#' \code{FUN(x)}.
#' @param x A fitted merMod object from \code{\link[lme4]{lmer}}.
#' @param FUN A function taking a fitted merMod object as input and returning
#'   the statistic of interest.
#' @param index An integer stating the position of the statistic in the output of
#'   \code{FUN(x)}.
#' @return A numeric vector with length equals to number of clusters of
#'   \code{x} containing the weighted influence value of each cluster.
#' @references Van der Leeden, R., Meijer, E., & Busing, F. M. T. A. (2008).
#'   Resampling multilevel models. In J. de Leeuw & E. Meijer (Eds.), Handbook
#'   of multilevel Analysis (pp. 401–433). New York, NY: Springer.
#' @export
#' @examples
#' library(lme4)
#' fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
#' # Define function for intraclass correlation
#' icc <- function(x) 1 / (1 + 1 / getME(x, "theta")^2)
#' empinf_mer(fm01ML, icc)
#' empinf_mer(fm01ML, fixef)
empinf_mer <- function(x, FUN, index = 1) {
  if (length(x@flist) > 1) {
    stop("currently can only compute influence values with two levels")
  }
  J <- lme4::ngrps(x)[[1]]
  th_noj <- rep(NA, J)
  n <- nobs(x)
  gp <- x@flist[[1]]
  gp_lv <- levels(gp)
  nj <- as.vector(unname(table(gp, dnn = NULL)))
  hj <- n / nj
  formula_x <- formula(x)
  org_data <- x@frame
  th_n <- FUN(x)[index]
  for (j in seq_along(th_noj)) {
    i <- which(gp == gp_lv[j])
    # m <- lmer(formula_x, data = org_data[-i, ])
    m_call <- update(x, formula_x, data = org_data[-i, ], evaluate = FALSE)
    m <- eval(m_call)
    th_noj[j] <- FUN(m)[index]
  }
  th_Jmj <- J * th_n - sum( (1 - 1 / hj) * th_noj)
  th_tilde <- hj * th_n - (hj - 1) * th_noj
  # (hj - 1) * (th_Jmj - th_tilde)
  (hj - 1) * (th_tilde - th_Jmj)
}

#' @rdname empinf_mer
#' @export
empinf_merm <- function(x, FUN) {
  if (length(x@flist) > 1) {
    stop("currently can only compute influence values with two levels")
  }
  J <- lme4::ngrps(x)[[1]]
  n <- nobs(x)
  gp <- x@flist[[1]]
  gp_lv <- levels(gp)
  nj <- as.vector(unname(table(gp, dnn = NULL)))
  hj <- n / nj
  formula_x <- formula(x)
  org_data <- x@frame
  th_n <- FUN(x)
  th_noj <- matrix(NA, nrow = J, ncol = length(th_n))
  for (j in seq_len(J)) {
    i <- which(gp == gp_lv[j])
    m_call <- update(x, formula_x, data = org_data[-i, ], evaluate = FALSE)
    m <- eval(m_call)
    th_noj[j, ] <- FUN(m)
  }
  th_Jmj <- J * th_n - colSums( (1 - 1 / hj) * th_noj)
  th_tilde <- tcrossprod(hj, th_n) - (hj - 1) * th_noj
  (hj - 1) * t(t(th_tilde) - th_Jmj)
}

empinf_mer_old <- function(x, FUN) {
  lj <- rep(NA, lme4::ngrps(x)[1])
  n <- nobs(x)
  gp <- x@flist[[1]]
  for (j in seq_along(lj)) {
    i <- which(gp == levels(gp)[j])
    nj <- length(i)
    m <- lmer(formula(x), data = x@frame[-i, ])
    lj[j] <- (1 - nj / n) * FUN(m)
  }
  FUN(x) - lj
}
