#' Completely Randomized Design
#'
#' Creates a CRD design object for random treatment assignment. Optionally
#' supports rerandomization for covariate balance.
#'
#' @param treated_ratio Numeric. Fraction of units assigned to treatment.
#'   Default 0.5.
#' @param covariate Numeric matrix. Covariates for balance checking (optional).
#' @param balance Logical. Use rerandomization for covariate balance. Default
#'   FALSE.
#' @param eps Numeric. Balance criterion threshold. Default 0.1.
#' @param max_iter Integer. Maximum rerandomization attempts. Default 1000.
#' @return A \code{crd} design object.
#' @export
crd <- function(treated_ratio = 0.5, covariate = NULL, balance = FALSE,
                eps = 0.1, max_iter = 1000) {
  obj <- list(
    params = treated_ratio,
    X = covariate,
    balance = balance,
    eps = eps,
    max_iter = max_iter
  )
  class(obj) <- c("crd", "design")
  obj
}

#' Bernoulli Randomized Design
#'
#' Creates a Bernoulli design where each unit is independently assigned to
#' treatment with probability \code{treated_prob}.
#'
#' @param treated_prob Numeric. Probability of treatment. Default 0.5.
#' @param covariate Numeric matrix. Covariates for balance checking (optional).
#' @param balance Logical. Use rerandomization. Default FALSE.
#' @param eps Numeric. Balance criterion threshold. Default 0.1.
#' @param max_iter Integer. Maximum rerandomization attempts. Default 1000.
#' @return A \code{bernoulli} design object.
#' @export
bernoulli <- function(treated_prob = 0.5, covariate = NULL, balance = FALSE,
                      eps = 0.1, max_iter = 1000) {
  obj <- list(
    params = treated_prob,
    X = covariate,
    balance = balance,
    eps = eps,
    max_iter = max_iter
  )
  class(obj) <- c("bernoulli", "design")
  obj
}

#' Infer Design Parameters from Observed Treatment
#'
#' @param design A design object.
#' @param Z Observed treatment vector.
#' @return Updated design object.
#' @export
get_params_via_obs <- function(design, Z) {
  UseMethod("get_params_via_obs")
}

#' @export
get_params_via_obs.crd <- function(design, Z) {
  design$params <- mean(Z)
  design
}

#' @export
get_params_via_obs.bernoulli <- function(design, Z) {
  design$params <- mean(Z)
  design
}

#' Draw Treatment Assignment
#'
#' @param design A design object.
#' @param n Number of units.
#' @return Numeric vector of treatment assignments (0/1).
#' @export
draw <- function(design, n) {
  UseMethod("draw")
}

#' @export
draw.crd <- function(design, n) {
  nt <- as.integer(design$params * n)
  template <- c(rep(1, nt), rep(0, n - nt))

  if (design$balance) {
    if (is.null(design$X)) {
      stop("covariate must be provided if balance is TRUE")
    }
    bal <- 1
    count <- 0
    Z <- template
    while (bal > design$eps && count < design$max_iter) {
      Z <- sample(template)
      bal <- get_balance(Z, design$X)
      count <- count + 1
    }
    if (bal > design$eps) {
      warning("Exceeded maximum iterations without achieving balance")
    }
    return(Z)
  }

  sample(template)
}

#' @export
draw.bernoulli <- function(design, n) {
  sample(c(0, 1), n, replace = TRUE)
}

#' Compute Rerandomization Balance Criterion
#'
#' Mahalanobis-type balance criterion from Morgan & Rubin (2015).
#' \deqn{m = \frac{n_t n_c}{n} (\bar{X}_t - \bar{X}_c)' \Sigma^{-1}
#'   (\bar{X}_t - \bar{X}_c)}
#'
#' @param Z Numeric vector. Treatment indicator.
#' @param X Numeric matrix. Covariates.
#' @return Numeric scalar. Balance criterion value.
#' @export
get_balance <- function(Z, X) {
  X <- as.matrix(X)
  n <- length(Z)
  idx_t <- Z == 1
  idx_c <- Z == 0
  nt <- sum(idx_t)
  nc <- sum(idx_c)

  m1 <- colMeans(X[idx_t, , drop = FALSE])
  m0 <- colMeans(X[idx_c, , drop = FALSE])
  Sigma <- stats::cov(X)

  diff <- m1 - m0
  m <- (nt * nc / n) * as.numeric(t(diff) %*% solve(Sigma) %*% diff)
  m
}
