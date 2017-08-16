#' Score Functions and Case-wise Derivatives
#'
#' @param x A fitted \code{merMod} object from \code{\link[lme4]{lmer}}.
#' @param level If \code{level = 1}, scores at level-1 are returned; if
#'   \code{level = 2}, which is the default, aggregated scores at the cluster-
#'   level are returned.
#' @export
scores_mer <- function(x, level = 2) {
  if (!lme4::isLMM(x)) {
    stop("currently only linear mixed model of class `merMod` is supported")
  }
  if (length(x@flist) > 1 & level == 2) {
    stop("set `level = 1` for models with multiple clustering levels")
  }
  # Extract required components
  PR <- x@pp
  X <- PR$X
  Lambdat <- PR$Lambdat
  Zt <- PR$Zt
  A <- Lambdat %*% Zt
  L <- PR$L()
  Lind <- PR$Lind
  RX <- PR$RX()
  n1 <- nobs(x)
  nfixed <- length(x@beta)
  Vinv <- Matrix::Diagonal(n = n1) -
    crossprod(Matrix::solve(L, Matrix::solve(L, A, system = "P"), system = "L"))
  VinvX <- Vinv %*% X
  eps <- x@resp$y - X %*% x@beta  # y - X %*% beta-hat
  sigma_x <- sigma(x)
  # Need REML
  if (lme4::isREML(x)) {
    B <- Vinv - crossprod(Matrix::solve(t(RX), t(VinvX)))
    df <- n1 - nfixed
  } else {
    B <- Vinv
    df <- n1
  }
  # score for beta
  score_beta <- VinvX * eps / sigma_x^2
  # d2_beta <- VinvX * X
  # score for sigma^2
  p1 <- (df / n1) / sigma_x^2
  p2 <- crossprod(Vinv, eps) * eps / sigma_x^4
  score_sigma2 <- (- p1 + p2) / 2
  # d2_sigma2 <- (p1 / 2 + p2) / sigma_x^2
  # score for variance components (not includings sigma^2)
  # score_vc <- d2_vc <- vector("list", length(x@theta))
  score_vc <- vector("list", length(x@theta))
  for (i in seq_along(score_vc)) {
    dLambdat <- Lambdat
    dLambdat@x[] <- Lind
    dG <- Matrix::forceSymmetric(dLambdat == i, uplo = "U")
    dV <- crossprod(Zt, dG %*% Zt)

    score_vc[[i]] <- - Matrix::diag(B %*% dV) / 2 +
      crossprod(Vinv %*% dV %*% Vinv, eps) * eps / 2 / sigma_x^2
    # d2_vc[[i]] <- Matrix::diag(B %*% dV %*% B %*% dV) / 2 +
    #   crossprod(Vinv %*% dV %*% Vinv %*% dV %*% Vinv, eps) * eps / sigma_x^2
  }
  score_out <- cbind(score_beta,
                     do.call(cbind, score_vc) / sigma_x^2, score_sigma2)
  # d2_out <- cbind(d2_beta, do.call(cbind, d2_vc), Residual = d2_sigma2)
  # Add names
  vc_names <- make_vcnames(lme4::VarCorr(x), sd_cor = FALSE)
  # colnames(score_out)[-(1:nfixed)] <-
  #   colnames(d2_out)[-(1:nfixed)] <- vc_names
  colnames(score_out)[- (1:nfixed)] <- vc_names
  if (level == 2) {
    Zt_clus <- as(x@flist[[1]], Class = "sparseMatrix")
    score_out <- Zt_clus %*% score_out
    # d2_out <- Zt_clus %*% d2_out
  }
  # Return a matrix
  # list(score = score_out, d2 = d2_out)
  score_out
}

# inf_jack <- function(x) {
#   eps <- 1e-6
#   w_org <- (x@resp$weights - eps)
#   gp <- m1@flist[[1]]
#   uniq_gp <- unique(gp)
#   lj <- numeric(length(uniq_gp))
#   for (i in seq_along(lj)) {
#     w <- w_org
#     this_gp <- seq_len(nobs(x))[gp == uniq_gp[i]]
#     w[this_gp] <- w[this_gp] + nobs(x) / length(this_gp) * eps
#     # up_x <- update(x, data = x@frame, weights = w)
#     up_x <- update_w(x, weights = w)
#     lj[i] <- (fixef(x)[2] - fixef(up_x)[2]) / eps
#   }
#   lj
# }
#
# update_w <- function (object, weights.)
# {
#   if (is.null(call <- getCall(object)))
#     stop("object should contain a 'call' component")
#   call$weights <- weights.
#   eval(call)
#   # if (length(extras) > 0) {
#   #   existing <- !is.na(match(names(extras), names(call)))
#   #   for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
#   #   if (any(!existing)) {
#   #     call <- c(as.list(call), extras[!existing])
#   #     call <- as.call(call)
#   #   }
#   # }
#   # if (evaluate) {
#   #   ff <- environment(formula(object))
#   #   pf <- parent.frame()
#   #   sf <- sys.frames()[[1]]
#   #   tryCatch(eval(call, envir = ff), error = function(e) {
#   #     tryCatch(eval(call, envir = sf), error = function(e) {
#   #       eval(call, pf)
#   #     })
#   #   })
#   # }
#   # else call
# }
