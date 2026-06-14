# ===========================================================================
# Tests for interference / clustered estimators
# Replicates Python tests from CausalModel/tests/test_est_via_aipw.py
# and CausalModel/tests/Interference.ipynb
# ===========================================================================

skip_slow_interference <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("CAUSALMODEL_SLOW_TESTS"), "1"),
    "set CAUSALMODEL_SLOW_TESTS=1 to run slow clustered interference simulations"
  )
}

# ---- Data generation: dimension checks (test_get_fixed_cluster) ----

test_that("generate_fixed_cluster produces correct dimensions", {
  set.seed(42)
  clusters <- 100
  group_struct <- c(2, 3)
  k <- 2
  data <- generate_fixed_cluster(clusters = clusters, group_struct = group_struct,
                                  k = k, tau = 1)
  n_units <- clusters * sum(group_struct)

  expect_equal(length(data$Y), n_units)
  expect_equal(length(data$Z), n_units)
  expect_equal(nrow(data$X), n_units)
  expect_equal(ncol(data$X), k)
  expect_equal(ncol(data$G), length(group_struct))
  expect_equal(length(data$cluster_labels), n_units)
  expect_equal(length(data$group_labels), n_units)
  expect_equal(length(data$ingroup_labels), n_units)
  expect_equal(ncol(data$Xc), 2 * k)
})

test_that("generate_fixed_cluster with 3 groups matches Python dimensions", {
  # Direct port of test_get_fixed_cluster from test_est_via_aipw.py
  set.seed(42)
  clusters <- 2
  group_struct <- c(3, 5, 7)
  k <- 11
  data <- generate_fixed_cluster(clusters = clusters, group_struct = group_struct,
                                  k = k, tau = 1,
                                  gamma = 0.01 * (0:(length(group_struct) - 1)))
  n_units <- clusters * sum(group_struct)

  expect_equal(length(data$Y), n_units)
  expect_equal(length(data$Z), n_units)
  expect_equal(nrow(data$X), n_units)
  expect_equal(ncol(data$X), k)
  expect_equal(ncol(data$G), length(group_struct))
  expect_equal(length(data$cluster_labels), n_units)
  expect_equal(length(data$group_labels), n_units)
  expect_equal(length(data$ingroup_labels), n_units)
  expect_equal(ncol(data$Xc), 2 * k)
})


# ---- Data generation: generate_clustered_data ----

test_that("generate_clustered_data works with multiple structures", {
  set.seed(42)
  clusters_list <- c(50, 50)
  group_struct_list <- list(c(2, 3), c(3, 4))
  data <- generate_clustered_data(
    clusters_list = clusters_list,
    group_struct_list = group_struct_list,
    tau = 1, gamma = c(0.1, 0.1)
  )
  n_units <- 50 * 5 + 50 * 7
  expect_equal(length(data$Y), n_units)
  expect_equal(length(data$Z), n_units)
  expect_equal(nrow(data$X), n_units)
})

test_that("generate_clustered_data with 3 structures matches Python", {
  # Port of test_get_clustered_data from test_est_via_aipw.py
  set.seed(42)
  clusters_list <- c(17, 19, 23)
  group_struct_list <- list(c(2, 3), c(5, 7), c(11, 13))
  data <- generate_clustered_data(
    clusters_list = clusters_list,
    group_struct_list = group_struct_list,
    tau = 1, gamma = c(0, 1)
  )
  k <- 2
  n_units <- sum(clusters_list * sapply(group_struct_list, sum))

  expect_equal(length(data$Y), n_units)
  expect_equal(length(data$Z), n_units)
  expect_equal(nrow(data$X), n_units)
  expect_equal(ncol(data$X), k)
  expect_equal(ncol(data$G), 2)
  expect_equal(length(data$cluster_labels), n_units)
  expect_equal(length(data$group_labels), n_units)
  expect_equal(length(data$ingroup_labels), n_units)
  expect_equal(ncol(data$Xc), 2 * k)
})


# ---- G encoding / decoding ----

test_that("encode/decode G are inverse operations", {
  group_struct <- c(2, 3, 4)
  for (g_enc in 0:5) {
    g <- decode_G(g_enc, group_struct)
    result <- encode_G_single(g, group_struct)
    expect_equal(result, g_enc)
  }
})

test_that("encode_G_mat matches encode_G_single for multiple rows", {
  group_struct <- c(3, 4)
  G <- rbind(c(0, 0), c(1, 2), c(2, 3), c(3, 4))
  encoded <- encode_G_mat(G, group_struct)
  expected <- sapply(1:nrow(G), function(i) encode_G_single(G[i, ], group_struct))
  expect_equal(encoded, expected)
})


# ---- Clustered IPW estimator ----

test_that("Clustered IPW estimator runs and returns valid structure", {
  skip_on_cran()
  set.seed(42)
  data <- generate_fixed_cluster(clusters = 500, group_struct = c(2, 3),
                                  k = 2, tau = 1, gamma = c(0.5, 0.1))
  cl <- clustered(data$Y, data$Z, data$X,
                   data$cluster_labels, data$group_labels,
                   data$ingroup_labels, n_matches = 20L)
  result <- est_via_ipw(cl)
  expect_true(is.list(result))
  expect_equal(length(result), 2)  # 2 groups
  expect_true(!is.null(result[[1]]$beta_g))
  expect_true(!is.null(result[[1]]$se))
})


# ---- Clustered AIPW estimator ----

test_that("Clustered AIPW estimator runs and returns valid structure", {
  skip_on_cran()
  set.seed(42)
  data <- generate_fixed_cluster(clusters = 500, group_struct = c(2, 3),
                                  k = 2, tau = 1, gamma = c(0.5, 0.1))
  cl <- clustered(data$Y, data$Z, data$X,
                   data$cluster_labels, data$group_labels,
                   data$ingroup_labels, n_matches = 20L)
  result <- est_via_aipw(cl)
  expect_true(is.list(result))
  expect_equal(length(result), 2)  # 2 groups
  expect_true(!is.null(result[[1]]$beta_g))
  expect_true(!is.null(result[[1]]$se))
  expect_true(!is.null(result[[2]]$beta_g))
  expect_true(!is.null(result[[2]]$se))
})


# ---- AIPW recovers beta(g): core correctness test ----
# Port of test_estimate_beta_g from test_est_via_aipw.py
# Verifies that clustered AIPW recovers the heterogeneous treatment effects
# beta(g) = tau + gamma * g across the neighborhood exposure grid.

test_that("Clustered AIPW recovers beta(g) within tolerance", {
  skip_slow_interference()
  set.seed(42)
  clusters_list <- c(2000, 3000, 4000)
  group_struct_list <- list(c(2, 2), c(3, 2), c(3, 3))
  tau <- 42
  gamma <- c(-50, 120)

  data <- generate_clustered_data(
    clusters_list = clusters_list,
    group_struct_list = group_struct_list,
    tau = tau, gamma = gamma
  )
  cl <- clustered(data$Y, data$Z, data$X,
                   data$cluster_labels, data$group_labels,
                   data$ingroup_labels)
  result <- est_via_aipw(cl)

  # Compute expected beta(g) on the maximal grid
  max_group_struct <- do.call(pmax, group_struct_list)  # c(3, 3)
  g1_vals <- 0:max_group_struct[1]
  g2_vals <- 0:max_group_struct[2]
  expected <- outer(g1_vals, g2_vals, function(g1, g2) tau + gamma[1] * g1 + gamma[2] * g2)

  # Group 1 (j=1): beta_g should match expected for g1 = 0..(max-1), all g2
  beta1 <- result[[1]]$beta_g
  # Reshape to matrix: (max_g1+1) x (max_g2+1)
  beta1_mat <- matrix(beta1, nrow = max_group_struct[1] + 1,
                       ncol = max_group_struct[2] + 1)
  # Compare interior values (last row excluded per Python test)
  for (g1 in 0:(max_group_struct[1] - 1)) {
    for (g2 in 0:max_group_struct[2]) {
      if (!is.nan(beta1_mat[g1 + 1, g2 + 1]) && beta1_mat[g1 + 1, g2 + 1] != 0) {
        expect_true(
          abs(beta1_mat[g1 + 1, g2 + 1] - expected[g1 + 1, g2 + 1]) /
            abs(expected[g1 + 1, g2 + 1]) < 0.10,
          label = sprintf("Group 1 beta(%d,%d): got %.2f, expected %.2f",
                          g1, g2, beta1_mat[g1 + 1, g2 + 1], expected[g1 + 1, g2 + 1])
        )
      }
    }
  }

  # Group 2 (j=2): beta_g should match expected for all g1, g2 = 0..(max-1)
  beta2 <- result[[2]]$beta_g
  beta2_mat <- matrix(beta2, nrow = max_group_struct[1] + 1,
                       ncol = max_group_struct[2] + 1)
  for (g1 in 0:max_group_struct[1]) {
    for (g2 in 0:(max_group_struct[2] - 1)) {
      if (!is.nan(beta2_mat[g1 + 1, g2 + 1]) && beta2_mat[g1 + 1, g2 + 1] != 0) {
        expect_true(
          abs(beta2_mat[g1 + 1, g2 + 1] - expected[g1 + 1, g2 + 1]) /
            abs(expected[g1 + 1, g2 + 1]) < 0.10,
          label = sprintf("Group 2 beta(%d,%d): got %.2f, expected %.2f",
                          g1, g2, beta2_mat[g1 + 1, g2 + 1], expected[g1 + 1, g2 + 1])
        )
      }
    }
  }
})


# ---- Asymptotic normality of clustered AIPW (simulation) ----
# Port of aipw_example.py / aipw_example.ipynb
# Verifies that (beta_hat(g) - beta(g)) / SE ~ N(0, 1)

test_that("Clustered AIPW estimates are approximately normal", {
  skip_slow_interference()
  set.seed(42)

  clusters_list <- c(1600, 2400, 3200)
  group_struct_list <- list(c(2, 2), c(3, 2), c(3, 3))
  tau <- 1.5
  gamma <- c(-0.5, 1.2)
  n_reps <- 200

  max_gs <- do.call(pmax, group_struct_list)
  expected <- outer(0:max_gs[1], 0:max_gs[2],
                     function(g1, g2) tau + gamma[1] * g1 + gamma[2] * g2)

  # Collect standardized estimates for group 1, g = (0, 0)
  z_stats <- numeric(n_reps)
  for (rep in seq_len(n_reps)) {
    data <- generate_clustered_data(
      clusters_list = clusters_list,
      group_struct_list = group_struct_list,
      tau = tau, gamma = gamma
    )
    cl <- clustered(data$Y, data$Z, data$X,
                     data$cluster_labels, data$group_labels,
                     data$ingroup_labels)
    result <- est_via_aipw(cl)
    beta_hat <- result[[1]]$beta_g[1]
    se_hat <- result[[1]]$se[1]
    if (se_hat > 0 && !is.nan(beta_hat)) {
      z_stats[rep] <- (beta_hat - expected[1, 1]) / se_hat
    } else {
      z_stats[rep] <- NA
    }
  }

  z_stats <- z_stats[!is.na(z_stats)]
  # Mean should be close to 0, variance close to 1
  expect_true(abs(mean(z_stats)) < 0.3,
              label = sprintf("mean(z) = %.3f, expected ~0", mean(z_stats)))
  expect_true(abs(var(z_stats) - 1) < 0.5,
              label = sprintf("var(z) = %.3f, expected ~1", var(z_stats)))
})


# ---- Coverage rate test ----
# From Interference.ipynb: checks 95% CI coverage is approximately 95%

test_that("Clustered AIPW achieves approximately 95% coverage", {
  skip_slow_interference()
  set.seed(42)

  clusters <- 500
  group_struct <- c(2, 3)
  tau <- 1
  gamma_vec <- c(0.5, 0.1)
  n_reps <- 200

  # Expected beta(g) for g=(0,0)
  expected_00 <- tau

  covers <- logical(n_reps)
  for (rep in seq_len(n_reps)) {
    data <- generate_fixed_cluster(clusters = clusters,
                                    group_struct = group_struct,
                                    k = 2, tau = tau, gamma = gamma_vec)
    cl <- clustered(data$Y, data$Z, data$X,
                     data$cluster_labels, data$group_labels,
                     data$ingroup_labels, n_matches = 20L)
    result <- est_via_aipw(cl)
    beta_hat <- result[[1]]$beta_g[1]
    se_hat <- result[[1]]$se[1]
    if (se_hat > 0 && !is.nan(beta_hat)) {
      covers[rep] <- abs(beta_hat - expected_00) < 1.96 * se_hat
    } else {
      covers[rep] <- NA
    }
  }
  coverage <- mean(covers, na.rm = TRUE)
  # Should be roughly 95% (allow 80-100% range for small simulation)
  expect_true(coverage > 0.80,
              label = sprintf("Coverage = %.1f%%, expected ~95%%", 100 * coverage))
})
