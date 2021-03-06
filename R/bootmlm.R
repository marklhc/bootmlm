#' bootmlm: A package for bootstrap resampling with multilevel data.
#'
#' Currently, the \code{\link[lme4]{bootMer}()} function in the \pkg{lme4}
#' package only implements the parametric bootstrap and a limited version of
#' semiparametric (or residual) bootstrap. The \pkg{bootmlm} package provides
#' the function \code{bootstrap_mer}, which performs various parametric,
#' residual, and case bootstrap resampling for fitted model objects with the
#' \pkg{lme4} package.
#'
#' @section Limitations:
#' \itemize{
#'   \item Currently only support multilevel models (a.k.a. linear mixed-effects
#'   models) fitted by \code{\link[lme4]{lmer}()} with the \pkg{lme4} package.
#'   Support for categorical outcome (i.e., generalized linear mixed-effects
#'   models fitted by \code{\link[lme4]{glmer}()}) and for models fitted with
#'   the \pkg{nlme} package will be added in the future.
#'   \item Random effect block bootstrap (\code{type = 'reb'}) and case
#'   bootstrap (\code{type = 'case'}) only support two-level models.
#'   \item Bias-corrected and accelerated bootstrap (using
#'   \code{\link{empinf_mer}()}) only supports two-level models.
#' }
#'
#' @docType package
#' @name bootmlm
NULL
