sample_re <- function(x, replace = TRUE, ...) {
  x[sample.int(length(x), replace = replace, ...)]
}

#' @importFrom Matrix crossprod tcrossprod t mean
#' @importFrom stats nobs resid runif sd sigma hatvalues model.matrix
.resid_resample <- function(x, nsim = 1, seed = NULL,
                            type = c("residual", "residual_cgr",
                                     "residual_trans"),
                            corrected = FALSE) {
  type <- match.arg(type)
  # residual bootstrap
  # if (!is.null(seed)) {
  if (!missing(seed)) {
    set.seed(seed)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1)
  }
  # RNGstate <- .Random.seed

  # Extract required quantities from the S4 object
  PR <- x@pp
  Zt <- PR$Zt
  X <- PR$X
  fixed <- X %*% x@beta

  # Js <- lme4::ngrps(x)
  # qs <- lengths(x@cnms)
  # nqs <- Js * qs
  # nqseq <- rep.int(seq_along(nqs), nqs)
  #
  # bstar <- get_reflate_b(x)
  # bstar_lst <- split(bstar, nqseq)
  # ml <- lapply(seq_along(bstar_lst),
  #              # easier to work with the transposed version
  #              function(i) matrix(bstar_lst[[i]], nrow = qs[i]))

  if (type == "residual_trans") {
    # Transform residuals
    r <- x@resp$y - fixed  # the variance of r may not be V
    V <- get_V(x)
    # var_r <- (V - X %*% Matrix::chol2inv(RX) %*% t(X)) * (V != 0)
    R <- Matrix::chol(V)
    if (corrected) {
      RX <- PR$RX()
      Vr <- V - X %*% Matrix::chol2inv(RX) %*% t(X)
      Vr[V == 0] <- 0
      Rr <- try(Matrix::chol(Vr), silent = TRUE)
      if (inherits(Rr, "try-error")) {
        Zeta <- get_zeta_eigen(r, Vr)
      } else {
        Zeta <- get_zeta(r, Rr)
      }
    } else {
      Zeta <- get_zeta(r, R)  # can use a corrected version
    }

    ss <- replicate(nsim, {
      as.vector(fixed +
                  crossprod(R, Zeta[sample.int(length(Zeta), replace = TRUE)]))
    },
    simplify = FALSE)
  } else if (type %in% c("residual", "residual_cgr")) {
    if (type == "residual") {
      ml <- get_reflate_b(x)
      estar <- get_reflate_e(x)
    } else if (type == "residual_cgr") {
      ml <- get_reflate_b_cgr(x)
      estar <- get_reflate_e_cgr(x)
    }
    ss <- replicate(nsim, {
      # faster!, and easier to read
      bstar_new <- unlist(lapply(ml, function(m) {
        m[ , sample.int(ncol(m), replace = TRUE)]
      }))

      as.vector(fixed + crossprod(Zt, bstar_new) +
                  estar[sample.int(length(estar), replace = TRUE)])
    },
    simplify = FALSE)
  }
  ss
}

.reb_resample <- function(x, nsim = 1, seed = NULL, scale = FALSE) {
  if (length(x@flist) > 1) {
    stop("currently REB bootstrap only support one level of clustering")
  }
  if (!missing(seed)) {
    set.seed(seed)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1)
  }

  PR <- x@pp
  Zt <- PR$Zt
  X <- PR$X
  fixed <- X %*% x@beta

  r <- x@resp$y - fixed
  rr <- get_reb_resid(x, Zt, r, scale = scale)
  ml <- rr$ml
  el <- rr$el

  ss <- replicate(nsim, {
    b_new <- unlist(lapply(ml, function(m)
      m[ , sample.int(ncol(m), replace = TRUE)]))
    e_new <- unlist(lapply(el, function(e)
      e[sample.int(length(e), replace = TRUE)]))

    as.vector(fixed + crossprod(Zt, b_new) + e_new)
  },
  simplify = FALSE)
  ss
}

case_newsample1 <- function(data, N, group, uniq_gp, gp_length, fname) {
  new_index2 <- c(sample_re(uniq_gp))
  new_index1 <- lapply(new_index2, function(i) seq_len(N)[group == i])
  group_length <- gp_length[new_index2]
  new_group <- rep(seq_along(new_index2), group_length)
  new_data <- data[unlist(new_index1), , drop = FALSE]
  new_data[fname] <- new_group
  rownames(new_data) <- NULL
  new_data
}

case_newsample2 <- function(data, N, group, uniq_gp, gp_length, fname) {
  new_index2 <- c(sample_re(uniq_gp, replace = TRUE))
  new_index1 <- lapply(new_index2,
                       function(i) sample_re(seq_len(N)[group == i],
                                             replace = TRUE))
  group_length <- gp_length[new_index2]
  new_group <- rep(seq_along(new_index2), group_length)
  new_data <- data[unlist(new_index1), , drop = FALSE]
  new_data[fname] <- new_group
  rownames(new_data) <- NULL
  new_data
}

.case_resample <- function(x, nsim = 1, seed = NULL,
                           lv1_resample = FALSE) {
  # if (!is.null(seed)) {
  if (length(x@flist) > 1) {
    stop("currently case bootstrap only support one level of clustering")
  }
  if (!missing(seed)) {
    set.seed(seed)
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    runif(1)
  }
  group <- as.numeric(x@flist[[1]])
  uniq_gp <- unique(group)
  gp_length <- unname(table(group))
  N <- nobs(x)
  org_data <- x@frame
  fname <- names(x@flist[1])

  if (!lv1_resample) {
    resample_fun <- case_newsample1
  } else {
    resample_fun <- case_newsample2
  }
  ss <- replicate(nsim,
                  resample_fun(org_data, N, group, uniq_gp, gp_length, fname),
                  simplify = FALSE)
  ss
}

# pupcross <- haven::read_sas(
#   "https://stats.idre.ucla.edu/wp-content/uploads/2016/02/pupcross.sas7bdat")
#
# library(lme4)
# x <- lmer(ACHIEV ~ PUPSEX + PUPSES + (PUPSES | PSCHOOL) + (1 | SSCHOOL),
#           data = pupcross)
# .resid_resample(x, 1000)
# .resid_cgr_resample(x, 1000)
# .resid_trans_resample(x, 1000)
