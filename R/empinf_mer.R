#' Empirical Influence Values for Two-Level Mixed Models
#'
#' This function calculates the empirical influence values for a statistic in
#'   a given fitted data object using the delete-\eqn{m_j} jackknife.
#'
#' \code{empinf_mer} computes non-parameteric influence function of models
#' fitted using \code{\link[lme4]{lmer}} by deleting one cluster at a time. See
#' van der Leeden, Meijer, and Busing (2008, pp. 420--422) for more information.
#' @param x A fitted merMod object from lmer.
#' @param FUN A function taking a fitted merMod object as input and returning
#'   the statistic of interest, which must be a (possibly named) numeric vector
#'   of length 1.
#' @param index An integer stating the position of the statistic in the output of
#'   \code{FUN(x)}.
#' @return A numeric vector with length equals to number of clusters of
#'   \code{x} containing the weighted influence value of each cluster.
#' @export
#' @examples
#' library(lme4)
#' fm01ML <- lmer(Yield ~ (1 | Batch), Dyestuff, REML = FALSE)
#' # Define function for intraclass correlation
#' icc <- function(x) 1 / (1 + 1 / getME(x, "theta")^2)
#' empinf_mer(fm01ML, icc)
empinf_mer <- function(x, FUN, index = 1) {
  if (length(x@flist) > 1) {
    stop("currently can only compute influence values with two levels")
  }
  J <- lme4::ngrps(x)[[1]]
  th_noj <- rep(NA, J)
  n <- nobs(x)
  gp <- x@flist[[1]]
  nj <- unname(table(gp))
  hj <- n / nj
  formula_x <- formula(x)
  org_data <- x@frame
  th_n <- FUN(x)[index]
  for (j in seq_along(th_noj)) {
    i <- which(gp == levels(gp)[j])
    m <- lmer(formula_x, data = org_data[-i, ])
    th_noj[j] <- FUN(m)[index]
  }
  th_Jmj <- J * th_n - sum((1 - 1 / hj) * th_noj)
  th_tilde <- hj * th_n - (hj - 1) * th_noj
  (hj - 1) * (th_Jmj - th_tilde)
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
  # lj <- lj - FUN(x)
  # sum(lj^3) / sum(lj^2)^1.5 / 6
  lj - FUN(x)
}

# Example (Not run):
# fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
# ff <- function(x) x@theta[1] / sigma(x)
# boo_fm1 <- bootMer(fm1, FUN = ff, nsim = 999)
# boot.ci(boo_fm1, type = c("norm", "basic", "perc"))
# boot.ci(boo_fm1, type = c("bca"), L = empinf.merMod(fm1, ff))
