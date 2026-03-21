#' Create a Causal Model Result
#'
#' Stores the output of a causal inference estimator: point estimate, standard
#' error, z-statistic, p-value, and 95% confidence interval.
#'
#' @param ate Numeric. Average treatment effect estimate.
#' @param se Numeric. Standard error of the ATE.
#' @param alpha Numeric. Significance level for CI (default 0.05).
#' @return A \code{cm_result} object (S3 list).
#' @export
cm_result <- function(ate, se, alpha = 0.05) {
  z_val <- ate / se
  p_value <- 2 * (1 - stats::pnorm(abs(z_val)))
  crit <- stats::qnorm(1 - alpha / 2)
  ci <- c(ate - crit * se, ate + crit * se)

  structure(
    list(
      ate = ate,
      se = se,
      z = z_val,
      p_value = p_value,
      ci = ci
    ),
    class = "cm_result"
  )
}

#' @export
print.cm_result <- function(x, ...) {
  cat("Causal Model Estimation Result\n")
  cat("------------------------------\n")
  cat(sprintf("  ATE:      %.6f\n", x$ate))
  cat(sprintf("  SE:       %.6f\n", x$se))
  cat(sprintf("  z:        %.4f\n", x$z))
  cat(sprintf("  p-value:  %.4f\n", x$p_value))
  cat(sprintf("  95%% CI:   [%.6f, %.6f]\n", x$ci[1], x$ci[2]))
  invisible(x)
}

#' @export
summary.cm_result <- function(object, ...) {
  data.frame(
    ATE = object$ate,
    SE = object$se,
    z = object$z,
    p_value = object$p_value,
    CI_lower = object$ci[1],
    CI_upper = object$ci[2]
  )
}
