#' Score Functions and Case-wise Derivatives
#'
#' @export
scores_mer <- function(x) {
  PR <- x@pp
  X <- PR$X
  Lambdat <- PR$Lambdat
  Zt <- PR$Zt
  A <- Lambdat %*% Zt
  L <- PR$L()
  RX <- PR$RX()
  n1 <- nobs(x)
  Vinv <- Matrix::Diagonal(n = n1) -
    crossprod(Matrix::solve(L, Matrix::solve(L, A, system = "P"), system = "L"))
  VinvX <- Vinv %*% X
  # VinvX <- X - crossprod(Matrix::solve(L, Matrix::solve(L, A, system = "P"),
  #                                      system = "L")) %*% X
  eps <- x@resp$y - X %*% x@beta
  sigma_x <- sigma(x)
  # Need REML
  if (lme4::isREML(x)) {
    B <- Vinv - crossprod(Matrix::solve(t(RX), t(VinvX)))
    df <- n1 - length(x@beta)
  } else {
    B <- Vinv
    df <- n1
  }

  p1 <- (df / n1) / sigma_x^2
  p2 <- crossprod(Vinv, eps) * eps / sigma_x^4
  score_sigma2 <- (- p1 + p2) / 2
  d2_sigma2 <- (p1 / 2 + p2) / sigma_x^2

  score_vc <- d2_vc <- vector("list", length(x@theta))
  for (i in seq_along(score_vc)) {
    dLambdat <- Lambdat
    dLambdat@x[] <- PR$Lind
    dG <- forceSymmetric(dLambdat == i, uplo = "U")
    dV <- crossprod(Zt, dG %*% Zt)

    score_vc[[i]] <- - Matrix::diag(B %*% dV) / 2 +
      crossprod(Vinv %*% dV %*% Vinv, eps) * eps / 2 / sigma_x^2
    d2_vc[[i]] <- Matrix::diag(B %*% dV %*% B %*% dV) / 2 +
      crossprod(Vinv %*% dV %*% Vinv %*% dV %*% Vinv, eps) * eps / sigma_x^2
  }

  score_beta <- VinvX * eps
  d2_beta <- VinvX * X

  list(score = cbind(score_beta, do.call(cbind, score_vc), score_sigma2),
       d2 = cbind(d2_beta, do.call(cbind, d2_vc), d2_sigma2))
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
#     lj[i] <- (fixef(up_x)[2] - fixef(x)[2]) / eps
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
