% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profile.R
\name{prof_ci_icc}
\alias{prof_ci_icc}
\title{Profile Likelihood Confidence Interval for Intraclass Correlation}
\usage{
prof_ci_icc(x, level = 0.95, eps_max = 1e+06)
}
\arguments{
\item{x}{A fitted merMod object from \code{\link[lme4]{lmer}}.}

\item{level}{Confidence level between 0 and 1. Default is .95.}

\item{eps_max}{The maximum value that the upper limit can be. Default is
1e6.}
}
\value{
A named numeric vector with two values showing the lower and upper
  limit of the intraclass correlation.
}
\description{
Compute confidence intervals for the intraclass correlation of a model
fit of class \code{\link[lme4]{merMod-class}}.
}
\details{
This function uses the \code{\link[stats]{uniroot}} function and determine
the lower and upper limit by evaluating the profile deviance (for ML) or the
-2 profile REML criterion. It works by obtaining the interval for
\eqn{\theta = \tau / \sigma} and transforming the two limits.

The resulting interval is bounded by zero for its lower limit.
}
