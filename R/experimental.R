#' Create an Experimental Study Object
#'
#' Constructs an experimental study data object for estimating causal effects
#' from randomized experiments.
#'
#' @param Y Numeric vector. Outcome variable.
#' @param Z Numeric vector. Binary treatment indicator (0/1).
#' @param X Numeric matrix, vector, or NULL. Covariates (optional).
#' @param design A design object (e.g., from \code{crd()} or
#'   \code{bernoulli()}). Default is CRD.
#' @return An \code{experimental} object (S3 list).
#' @export
experimental <- function(Y, Z, X = NULL, design = NULL) {
  if (is.null(X)) {
    X <- matrix(rep(1, length(Y)), ncol = 1)
  }
  data <- po_data(Y, Z, X)

  if (is.null(design)) design <- crd()
  design <- get_params_via_obs(design, Z)

  obj <- list(
    data = data,
    design = design,
    stats = NULL,
    cal_stats = NULL
  )
  class(obj) <- "experimental"
  obj
}

#' @rdname estimate
#' @export
estimate.experimental <- function(obj, ...) {
  est_via_dm(obj, ...)
}

#' Difference-in-Means Estimator
#'
#' Simple difference in sample means between treated and control groups.
#'
#' @param obj An \code{experimental} object.
#' @param ... Additional arguments (currently unused).
#' @return A \code{cm_result} object.
#' @export
est_via_dm <- function(obj, ...) {
  UseMethod("est_via_dm")
}

#' @rdname est_via_dm
#' @export
est_via_dm.experimental <- function(obj, ...) {
  d <- obj$data
  dm <- cal_dm(d$Z, d$Y)
  cm_result(dm$ate, dm$se)
}

#' Stratified Estimator
#'
#' Estimates ATE as a weighted average of within-stratum difference-in-means.
#'
#' @param obj An \code{experimental} object.
#' @param ... Additional arguments passed to methods.
#' @return A \code{cm_result} object.
#' @export
est_via_strata <- function(obj, ...) {
  UseMethod("est_via_strata")
}

#' @rdname est_via_strata
#' @param strata Integer or factor vector of stratum labels (same length as Y).
#' @export
est_via_strata.experimental <- function(obj, strata, ...) {
  d <- obj$data
  if (length(strata) != d$n) {
    stop("strata must have the same length as Y")
  }

  strata_levels <- unique(strata)
  ate_parts <- numeric(length(strata_levels))
  se2_parts <- numeric(length(strata_levels))

  for (i in seq_along(strata_levels)) {
    lev <- strata_levels[i]
    mask <- strata == lev
    w <- mean(mask)
    dm <- cal_dm(d$Z[mask], d$Y[mask])
    ate_parts[i] <- dm$ate * w
    se2_parts[i] <- (dm$se^2) * (w^2)
  }

  ate <- sum(ate_parts)
  se <- sqrt(sum(se2_parts))
  cm_result(ate, se)
}

#' ANCOVA Estimator
#'
#' Fisher's ANCOVA with treatment-covariate interactions and HC0 robust
#' standard errors. Regresses Y ~ 1 + Z + X + Z*X.
#'
#' @param obj An \code{experimental} object.
#' @param ... Additional arguments (currently unused).
#' @return A \code{cm_result} object.
#' @export
est_via_ancova <- function(obj, ...) {
  UseMethod("est_via_ancova")
}

#' @rdname est_via_ancova
#' @export
est_via_ancova.experimental <- function(obj, ...) {
  d <- obj$data
  Z_col <- matrix(d$Z, ncol = 1)
  regressor <- cbind(1, Z_col, d$X, d$X * d$Z)

  XtX_inv <- solve(crossprod(regressor))
  beta <- XtX_inv %*% crossprod(regressor, d$Y)
  resid <- d$Y - regressor %*% beta

  meat <- crossprod(regressor * as.numeric(resid^2), regressor)
  vcov_hc0 <- XtX_inv %*% meat %*% XtX_inv

  ate <- beta[2]
  se <- sqrt(vcov_hc0[2, 2])
  cm_result(ate, se)
}

#' Fisher Randomization Test
#'
#' Computes a p-value by comparing the observed test statistic (DM) to
#' the distribution under random treatment reassignment.
#'
#' @param obj An \code{experimental} object.
#' @param ... Additional arguments passed to methods.
#' @return Numeric p-value.
#' @export
test_via_fisher <- function(obj, ...) {
  UseMethod("test_via_fisher")
}

#' @rdname test_via_fisher
#' @param n_draws Number of random draws. Default 1000.
#' @export
test_via_fisher.experimental <- function(obj, n_draws = 1000, ...) {
  d <- obj$data
  obs_stat <- mean(d$Y[d$Z == 1]) - mean(d$Y[d$Z == 0])

  cal_stat <- function(Z_s) {
    mean(d$Y[Z_s == 1]) - mean(d$Y[Z_s == 0])
  }

  T_s <- numeric(n_draws)
  for (i in seq_len(n_draws)) {
    Z_s <- draw(obj$design, d$n)
    T_s[i] <- cal_stat(Z_s)
  }

  pval <- min(mean(T_s > obs_stat), mean(T_s < obs_stat))
  pval
}


# ---- Internal helpers ----

#' @keywords internal
cal_dm <- function(Z, Y) {
  Y1 <- Y[Z == 1]
  Y0 <- Y[Z == 0]
  ate <- mean(Y1) - mean(Y0)
  v1 <- mean((Y1 - mean(Y1))^2)
  v0 <- mean((Y0 - mean(Y0))^2)
  se <- sqrt(v1 / sum(Z == 1) + v0 / sum(Z == 0))
  list(ate = ate, se = se)
}
