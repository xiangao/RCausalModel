test_that("Difference-in-means estimator works", {
  set.seed(42)
  N <- 2000
  tau <- 5
  Y0 <- rnorm(N)
  Z <- rep(c(0, 1), each = N / 2)
  Y <- Y0 + tau * Z
  exp_obj <- experimental(Y, Z)
  result <- est_via_dm(exp_obj)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - tau) < 1)
})

test_that("Stratified estimator works", {
  set.seed(42)
  N <- 2000
  tau <- 5
  strata <- rep(c(0, 1), each = N / 2)
  Y0 <- rnorm(N) + 2 * strata
  Z <- rep(c(0, 1, 0, 1), each = N / 4)
  Y <- Y0 + tau * Z
  exp_obj <- experimental(Y, Z)
  result <- est_via_strata(exp_obj, strata = strata)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - tau) < 1.5)
})

test_that("ANCOVA estimator works", {
  set.seed(42)
  N <- 2000
  tau <- 5
  X <- matrix(rnorm(N * 2), ncol = 2)
  Y0 <- X %*% c(1, 2) + rnorm(N)
  Z <- sample(c(0, 1), N, replace = TRUE)
  Y <- as.numeric(Y0) + tau * Z
  exp_obj <- experimental(Y, Z, X)
  result <- est_via_ancova(exp_obj)
  expect_s3_class(result, "cm_result")
  expect_true(abs(result$ate - tau) < 1.5)
})

test_that("Fisher randomization test returns valid p-value", {
  set.seed(42)
  N <- 200
  tau <- 3
  Y0 <- rnorm(N)
  Z <- rep(c(0, 1), each = N / 2)
  Y <- Y0 + tau * Z
  exp_obj <- experimental(Y, Z)
  pval <- test_via_fisher(exp_obj, n_draws = 500)
  expect_true(is.numeric(pval))
  expect_true(pval >= 0 && pval <= 1)
  # With tau=3 and N=200, should reject
  expect_true(pval < 0.05)
})

test_that("CRD design draws correct number of treated", {
  design <- crd(treated_ratio = 0.4)
  Z <- draw(design, 100)
  expect_equal(sum(Z), 40)
  expect_equal(length(Z), 100)
})

test_that("Bernoulli design works", {
  set.seed(42)
  design <- bernoulli(treated_prob = 0.5)
  Z <- draw(design, 1000)
  expect_equal(length(Z), 1000)
  expect_true(all(Z %in% c(0, 1)))
  # Roughly half should be treated
  expect_true(abs(mean(Z) - 0.5) < 0.1)
})

test_that("Balance criterion is numeric", {
  X <- matrix(rnorm(200), ncol = 2)
  Z <- sample(c(0, 1), 100, replace = TRUE)
  b <- get_balance(Z, X)
  expect_true(is.numeric(b))
  expect_true(b >= 0)
})
