#' Create an Observational Study Object
#'
#' Constructs an observational study data object for estimating causal effects
#' from non-experimental data.
#'
#' @param Y Numeric vector. Outcome variable.
#' @param Z Numeric vector. Treatment variable (binary or continuous).
#' @param X Numeric matrix or vector. Covariates.
#' @return An \code{observational} object (S3 list).
#' @export
#' @examples
#' data <- generate_data(N = 1000, k = 2, tau = 5)
#' obs <- observational(data$Y, data$Z, data$X)
#' est_via_aipw(obs)
observational <- function(Y, Z, X) {
  data <- po_data(Y, Z, X)
  obj <- list(
    data = data,
    propensity = NULL,
    treated_pred = NULL,
    control_pred = NULL,
    eps = 1e-4
  )
  class(obj) <- "observational"
  obj
}

#' Default Estimation
#'
#' Uses AIPW for binary treatment and DML for continuous treatment.
#'
#' @param obj An \code{observational} or \code{experimental} object.
#' @param ... Additional arguments passed to the selected estimator.
#' @return A \code{cm_result} object.
#' @export
estimate <- function(obj, ...) {
  UseMethod("estimate")
}

#' @rdname estimate
#' @export
estimate.observational <- function(obj, ...) {
  if (length(unique(obj$data$Z)) == 2) {
    est_via_aipw(obj, ...)
  } else {
    est_via_dml(obj, ...)
  }
}

#' OLS Estimator
#'
#' Estimates ATE via linear regression of Y on (Z, X) with HC0
#' heteroskedasticity-robust standard errors.
#'
#' @param obj An \code{observational} object.
#' @param ... Additional arguments (currently unused).
#' @return A \code{cm_result} object.
#' @export
est_via_ols <- function(obj, ...) {
  UseMethod("est_via_ols")
}

#' @rdname est_via_ols
#' @export
est_via_ols.observational <- function(obj, ...) {
  d <- obj$data
  regressor <- cbind(d$Z, d$X)

  XtX_inv <- solve(crossprod(regressor))
  beta <- XtX_inv %*% crossprod(regressor, d$Y)
  resid <- d$Y - regressor %*% beta

  # HC0 robust standard errors
  meat <- crossprod(regressor * as.numeric(resid^2), regressor)
  vcov_hc0 <- XtX_inv %*% meat %*% XtX_inv

  ate <- beta[1]
  se <- sqrt(vcov_hc0[1, 1])

  cm_result(ate, se)
}

#' IPW Estimator
#'
#' Inverse Probability Weighting estimator for ATE.
#'
#' @param obj An \code{observational} or \code{clustered} object.
#' @param ... Additional arguments passed to methods.
#' @return A \code{cm_result} object.
#' @export
est_via_ipw <- function(obj, ...) {
  UseMethod("est_via_ipw")
}

#' @rdname est_via_ipw
#' @param propensity_model A model wrapper with \code{fit} and
#'   \code{predict_proba}. Default uses logistic regression.
#' @param propensity Optional numeric vector of pre-estimated propensity scores.
#' @param normalize Logical. Normalize weights (Hajek estimator). Default TRUE.
#' @export
est_via_ipw.observational <- function(obj, propensity_model = NULL,
                                       propensity = NULL, normalize = TRUE, ...) {
  d <- obj$data

  if (!is.null(propensity)) {
    ps <- propensity
  } else {
    if (is.null(propensity_model)) propensity_model <- cm_logistic()
    fitted_model <- propensity_model$fit(d$X, d$Z)
    ps <- propensity_model$insample_proba(fitted_model)
  }

  ps <- fix_propensity(ps, obj$eps)

  w1 <- d$Z / ps
  w0 <- (1 - d$Z) / (1 - ps)

  if (normalize) {
    G <- w1 * d$Y / (sum(w1) / d$n) - w0 * d$Y / (sum(w0) / d$n)
  } else {
    G <- w1 * d$Y - w0 * d$Y
  }

  ate <- mean(G)
  se <- sqrt(stats::var(G) / (length(G) - 1))
  cm_result(ate, se)
}

#' AIPW (Doubly Robust) Estimator
#'
#' Augmented IPW estimator. Consistent if either the outcome model or the
#' propensity model is correctly specified.
#'
#' @param obj An \code{observational} or \code{clustered} object.
#' @param ... Additional arguments passed to methods.
#' @return A \code{cm_result} object.
#' @export
est_via_aipw <- function(obj, ...) {
  UseMethod("est_via_aipw")
}

#' @rdname est_via_aipw
#' @param outcome_model A model wrapper with \code{fit} and \code{predict}.
#'   Default uses OLS.
#' @param propensity_model A model wrapper with \code{fit} and
#'   \code{predict_proba}. Default uses logistic regression.
#' @param treated_pred Optional pre-computed treated outcome predictions.
#' @param control_pred Optional pre-computed control outcome predictions.
#' @param propensity Optional pre-computed propensity scores.
#' @export
est_via_aipw.observational <- function(obj, outcome_model = NULL,
                                        propensity_model = NULL,
                                        treated_pred = NULL,
                                        control_pred = NULL,
                                        propensity = NULL, ...) {
  d <- obj$data
  if (is.null(outcome_model)) outcome_model <- cm_ols()
  if (is.null(propensity_model)) propensity_model <- cm_logistic()

  if (!is.null(treated_pred)) {
    mu1 <- treated_pred
  } else {
    fitted_t <- outcome_model$fit(d$Xt, d$Yt)
    mu1 <- outcome_model$predict(fitted_t, d$X)
  }

  if (!is.null(control_pred)) {
    mu0 <- control_pred
  } else {
    fitted_c <- outcome_model$fit(d$Xc, d$Yc)
    mu0 <- outcome_model$predict(fitted_c, d$X)
  }

  if (!is.null(propensity)) {
    ps <- propensity
  } else {
    fitted_ps <- propensity_model$fit(d$X, d$Z)
    ps <- propensity_model$insample_proba(fitted_ps)
  }

  ps <- fix_propensity(ps, obj$eps)

  G <- (mu1 - mu0) +
    d$Z * (d$Y - mu1) / ps -
    (1 - d$Z) * (d$Y - mu0) / (1 - ps)

  ate <- mean(G)
  se <- sqrt(stats::var(G) / (length(G) - 1))
  cm_result(ate, se)
}

#' Matching Estimator
#'
#' k-nearest-neighbor matching estimator for ATE, with optional bias
#' adjustment.
#'
#' @param obj An \code{observational} object.
#' @param ... Additional arguments passed to methods.
#' @return A \code{cm_result} object.
#' @export
est_via_matching <- function(obj, ...) {
  UseMethod("est_via_matching")
}

#' @rdname est_via_matching
#' @param num_matches Number of matches per unit (M). Default 1.
#' @param num_matches_for_var Number of matches for variance estimation (J).
#'   Default equals \code{num_matches}.
#' @param bias_adj Logical. Apply bias correction using OLS. Default FALSE.
#' @export
est_via_matching.observational <- function(obj, num_matches = 1,
                                            num_matches_for_var = NULL,
                                            bias_adj = FALSE, ...) {
  d <- obj$data
  M <- num_matches
  J <- if (is.null(num_matches_for_var)) M else num_matches_for_var

  sd_Xt <- apply(d$Xt, 2, stats::sd)
  sd_Xc <- apply(d$Xc, 2, stats::sd)
  sd_Xt[sd_Xt == 0] <- 1
  sd_Xc[sd_Xc == 0] <- 1
  Xt_scaled <- sweep(d$Xt, 2, sd_Xt, "/")
  Xc_scaled <- sweep(d$Xc, 2, sd_Xc, "/")

  match_for_t <- mat_match_mat(Xt_scaled, Xc_scaled, M)
  match_for_c <- mat_match_mat(Xc_scaled, Xt_scaled, M)

  Yhat_t <- rowMeans(matrix(d$Yc[match_for_t], nrow = d$nt, ncol = M))
  Yhat_c <- rowMeans(matrix(d$Yt[match_for_c], nrow = d$nc, ncol = M))
  ITT_t <- d$Yt - Yhat_t
  ITT_c <- Yhat_c - d$Yc

  att <- mean(ITT_t)
  atc <- mean(ITT_c)
  ate <- (d$nc / d$n) * atc + (d$nt / d$n) * att

  if (bias_adj) {
    ols <- cm_ols()
    mu0_model <- ols$fit(d$Xc, d$Yc)
    mu1_model <- ols$fit(d$Xt, d$Yt)
    mu0_t <- ols$predict(mu0_model, d$Xt)
    mu0_c <- ols$predict(mu0_model, d$Xc)
    mu1_t <- ols$predict(mu1_model, d$Xt)
    mu1_c <- ols$predict(mu1_model, d$Xc)
    match_for_0t <- rowMeans(matrix(mu0_c[match_for_t], nrow = d$nt, ncol = M))
    match_for_1c <- rowMeans(matrix(mu1_t[match_for_c], nrow = d$nc, ncol = M))
    BM <- sum(mu0_t - match_for_0t) / d$n - sum(mu1_c - match_for_1c) / d$n
    ate <- ate - BM
  }

  Yhat1 <- c(d$Yt, Yhat_c)
  Yhat0 <- c(Yhat_t, d$Yc)

  Km <- numeric(d$n)
  for (i in seq_len(nrow(match_for_c))) {
    for (j in seq_len(ncol(match_for_c))) {
      Km[match_for_c[i, j]] <- Km[match_for_c[i, j]] + 1
    }
  }
  for (i in seq_len(nrow(match_for_t))) {
    for (j in seq_len(ncol(match_for_t))) {
      Km[match_for_t[i, j] + d$nt] <- Km[match_for_t[i, j] + d$nt] + 1
    }
  }

  match_tt <- mat_match_mat(Xt_scaled, Xt_scaled, J + 1)
  match_cc <- mat_match_mat(Xc_scaled, Xc_scaled, J + 1)

  Yhat_tt <- rowMeans(matrix(d$Yt[match_tt], nrow = d$nt, ncol = J + 1))
  Yhat_cc <- rowMeans(matrix(d$Yc[match_cc], nrow = d$nc, ncol = J + 1))

  Y_all <- c(d$Yt, d$Yc)
  Y_close <- c(Yhat_tt, Yhat_cc)

  sigmaXW <- (J + 1) / J * (Y_all - Y_close)^2

  V1 <- (Yhat1 - Yhat0 - ate)^2
  V2 <- ((Km / M)^2 + (2 * M - 1) / M * Km / M) * sigmaXW
  V <- mean(V1 + V2)

  se <- sqrt(V / d$n)
  cm_result(ate, se)
}

#' Double/Debiased Machine Learning Estimator
#'
#' DML estimator for continuous treatments using cross-fitting.
#'
#' @param obj An \code{observational} object.
#' @param ... Additional arguments passed to methods.
#' @return A \code{cm_result} object.
#' @export
est_via_dml <- function(obj, ...) {
  UseMethod("est_via_dml")
}

#' @rdname est_via_dml
#' @param outcome_model Model wrapper for Y ~ X. Default uses OLS.
#' @param treatment_model Model wrapper for Z ~ X. Default uses OLS.
#' @param k_folds Number of cross-fitting folds. Default 2.
#' @export
est_via_dml.observational <- function(obj, outcome_model = NULL,
                                       treatment_model = NULL,
                                       k_folds = 2, ...) {
  d <- obj$data
  if (is.null(outcome_model)) outcome_model <- cm_ols()
  if (is.null(treatment_model)) treatment_model <- cm_ols()

  n <- d$n
  fold_ids <- sample(rep(seq_len(k_folds), length.out = n))

  thetas <- numeric(k_folds)
  phi2s <- numeric(k_folds)
  Js <- numeric(k_folds)

  for (k in seq_len(k_folds)) {
    test_idx <- which(fold_ids == k)
    train_idx <- which(fold_ids != k)

    om <- outcome_model$fit(d$X[train_idx, , drop = FALSE], d$Y[train_idx])
    U <- d$Y[test_idx] - outcome_model$predict(om, d$X[test_idx, , drop = FALSE])

    tm <- treatment_model$fit(d$X[train_idx, , drop = FALSE], d$Z[train_idx])
    V <- d$Z[test_idx] - treatment_model$predict(tm, d$X[test_idx, , drop = FALSE])

    theta <- sum(V * U) / sum(V * V)
    thetas[k] <- theta

    phi2s[k] <- mean(V^2 * (U - V * theta)^2)
    Js[k] <- mean(V^2)
  }

  ate <- mean(thetas)
  se <- sqrt(mean(phi2s) / (mean(Js)^2)) / sqrt(n)
  cm_result(ate, se)
}


# ---- Internal helpers ----

#' @keywords internal
fix_propensity <- function(ps, eps = 1e-4) {
  bad <- sum(ps * (1 - ps) == 0)
  if (bad > 0) {
    ps[ps == 0] <- eps
    ps[ps == 1] <- 1 - eps
    warning(sprintf("Propensity scores had %d values of 0 or 1, trimmed to [%g, %g].",
                    bad, eps, 1 - eps))
  }
  ps
}

#' @keywords internal
mat_match_mat <- function(X, Y, M) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n_Y <- nrow(Y)
  k <- min(M, n_Y)
  nn <- RANN::nn2(data = Y, query = X, k = k)
  idx <- nn$nn.idx

  if (M > n_Y) {
    extra <- M - n_Y
    pad <- matrix(sample(seq_len(n_Y), nrow(X) * extra, replace = TRUE),
                  nrow = nrow(X), ncol = extra)
    idx <- cbind(idx, pad)
  }

  idx
}
