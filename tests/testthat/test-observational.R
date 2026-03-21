test_that("OLS estimator recovers ATE", {
  set.seed(42)
  data <- generate_data(N = 5000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_ols(obs)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 2)
})

test_that("IPW estimator recovers ATE", {
  set.seed(42)
  data <- generate_data(N = 5000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_ipw(obs)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 2)
})

test_that("AIPW estimator recovers ATE", {
  set.seed(42)
  data <- generate_data(N = 5000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_aipw(obs)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 2)
})

test_that("Matching estimator recovers ATE", {
  set.seed(42)
  data <- generate_data(N = 5000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_matching(obs, num_matches = 3)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 3)
})

test_that("Matching with bias adjustment works", {
  set.seed(42)
  data <- generate_data(N = 3000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_matching(obs, num_matches = 3, bias_adj = TRUE)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 3)
})

test_that("DML estimator recovers ATE for continuous treatment", {
  set.seed(42)
  data <- generate_data_continuous(N = 5000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_dml(obs)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 2)
})

test_that("Default estimate dispatches correctly", {
  set.seed(42)
  # Binary Z -> AIPW
  data <- generate_data(N = 2000, k = 2, tau = 5)
  obs <- observational(data$Y, data$Z, data$X)
  result <- estimate(obs)
  expect_s3_class(result, "cm_result")

  # Continuous Z -> DML
  data_c <- generate_data_continuous(N = 2000, k = 2, tau = 5)
  obs_c <- observational(data_c$Y, data_c$Z, data_c$X)
  result_c <- estimate(obs_c)
  expect_s3_class(result_c, "cm_result")
})

# ---- From observational.ipynb: higher-dimensional covariates ----

test_that("AIPW works with k=5 covariates", {
  set.seed(42)
  data <- generate_data(N = 5000, k = 5, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_aipw(obs)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 3)
})

test_that("Matching with num_matches=10 recovers ATE", {
  set.seed(42)
  data <- generate_data(N = 5000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_matching(obs, num_matches = 10)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 3)
})

# ---- From observational.ipynb: DML with explicit k_folds ----

test_that("DML with k_folds=5 recovers ATE", {
  set.seed(42)
  data <- generate_data_continuous(N = 5000, k = 2, tau = 10)
  obs <- observational(data$Y, data$Z, data$X)
  result <- est_via_dml(obs, k_folds = 5)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - 10) < 2)
})

# ---- From observational.ipynb: asymptotic normality simulation ----

test_that("AIPW standardized estimates are approximately normal", {
  skip_on_cran()
  skip("Slow simulation test — run manually with test_file()")
  set.seed(42)
  n_reps <- 500
  tau <- 10
  z_stats <- numeric(n_reps)

  for (rep in seq_len(n_reps)) {
    data <- generate_data(N = 1000, k = 2, tau = tau)
    obs <- observational(data$Y, data$Z, data$X)
    result <- est_via_aipw(obs)
    z_stats[rep] <- (result$ate - tau) / result$se
  }

  # (tau_hat - tau) / SE should be ~ N(0, 1)
  expect_true(abs(mean(z_stats)) < 0.15,
              label = sprintf("mean(z) = %.3f, expected ~0", mean(z_stats)))
  expect_true(abs(var(z_stats) - 1) < 0.3,
              label = sprintf("var(z) = %.3f, expected ~1", var(z_stats)))
})

test_that("DML standardized estimates are approximately normal", {
  skip_on_cran()
  skip("Slow simulation test — run manually with test_file()")
  set.seed(42)
  n_reps <- 500
  tau <- 10
  z_stats <- numeric(n_reps)

  for (rep in seq_len(n_reps)) {
    data <- generate_data_continuous(N = 1000, k = 2, tau = tau)
    obs <- observational(data$Y, data$Z, data$X)
    result <- est_via_dml(obs, k_folds = 2)
    z_stats[rep] <- (result$ate - tau) / result$se
  }

  expect_true(abs(mean(z_stats)) < 0.15,
              label = sprintf("mean(z) = %.3f, expected ~0", mean(z_stats)))
  expect_true(abs(var(z_stats) - 1) < 0.3,
              label = sprintf("var(z) = %.3f, expected ~1", var(z_stats)))
})

# ---- Utilities ----

test_that("cm_result print and summary work", {
  r <- cm_result(5.0, 0.5)
  expect_output(print(r), "ATE")
  s <- summary(r)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 1)
})
