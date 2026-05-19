# Getting Started with CausalModel

## Overview

CausalModel provides causal inference estimators for three settings:

1.  **Observational studies** — OLS, IPW, AIPW, matching, DML
2.  **Randomized experiments** — difference-in-means, stratified,
    ANCOVA, Fisher test
3.  **Network interference** — clustered IPW and AIPW

This vignette walks through basic usage of each.

## Observational Study

### Generate synthetic data

The
[`generate_data()`](https://xiangao.github.io/RCausalModel/reference/generate_data.md)
function creates (Y, Z, X) with a known treatment effect `tau`.
Treatment assignment depends on covariates through a logistic propensity
model, so naive comparisons are confounded.

``` r

data <- generate_data(N = 2000, k = 2, tau = 10)
obs <- observational(data$Y, data$Z, data$X)
```

### OLS

Linear regression of Y on (Z, X) with HC0 robust standard errors:

``` r

est_via_ols(obs)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      9.982300
#>   SE:       0.034704
#>   z:        287.6389
#>   p-value:  0.0000
#>   95% CI:   [9.914281, 10.050319]
```

### IPW (Inverse Probability Weighting)

Reweights outcomes by inverse propensity scores. Uses logistic
regression by default:

``` r

est_via_ipw(obs)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      9.954672
#>   SE:       0.282654
#>   z:        35.2185
#>   p-value:  0.0000
#>   95% CI:   [9.400679, 10.508664]
```

### AIPW (Doubly Robust)

Combines outcome modeling and propensity weighting — consistent if
*either* model is correctly specified:

``` r

est_via_aipw(obs)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      9.951175
#>   SE:       0.052355
#>   z:        190.0694
#>   p-value:  0.0000
#>   95% CI:   [9.848560, 10.053789]
```

### Matching

Nearest-neighbor matching with `M` matches per unit. Variance estimation
follows Abadie & Imbens (2006):

``` r

est_via_matching(obs, num_matches = 3)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      10.043697
#>   SE:       0.057497
#>   z:        174.6809
#>   p-value:  0.0000
#>   95% CI:   [9.931004, 10.156390]
```

Bias-adjusted matching corrects for residual imbalance in matched
covariates:

``` r

est_via_matching(obs, num_matches = 3, bias_adj = TRUE)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      9.974851
#>   SE:       0.057518
#>   z:        173.4214
#>   p-value:  0.0000
#>   95% CI:   [9.862117, 10.087584]
```

### DML (Double Machine Learning)

For continuous treatments, DML with cross-fitting estimates the linear
treatment effect. First, generate continuous treatment data:

``` r

data_c <- generate_data_continuous(N = 2000, k = 2, tau = 10)
obs_c <- observational(data_c$Y, data_c$Z, data_c$X)
est_via_dml(obs_c)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      9.986547
#>   SE:       0.021793
#>   z:        458.2498
#>   p-value:  0.0000
#>   95% CI:   [9.943834, 10.029260]
```

### Default estimator

The
[`estimate()`](https://xiangao.github.io/RCausalModel/reference/estimate.md)
generic dispatches to AIPW for binary treatment and DML for continuous:

``` r

estimate(obs)    # calls est_via_aipw
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      9.951175
#>   SE:       0.052355
#>   z:        190.0694
#>   p-value:  0.0000
#>   95% CI:   [9.848560, 10.053789]
estimate(obs_c)  # calls est_via_dml
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      9.987542
#>   SE:       0.021876
#>   z:        456.5507
#>   p-value:  0.0000
#>   95% CI:   [9.944666, 10.030419]
```

## Randomized Experiment

In a completely randomized design, treatment is independent of
covariates. We simulate this directly:

``` r

N <- 1000
X <- matrix(rnorm(N * 2), ncol = 2)
tau <- 5
Z <- sample(c(rep(1, N/2), rep(0, N/2)))
Y <- tau * Z + X %*% c(1, -1) + rnorm(N)
exp_obj <- experimental(Y, Z, X)
```

### Difference-in-Means

``` r

est_via_dm(exp_obj)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      4.993835
#>   SE:       0.107146
#>   z:        46.6076
#>   p-value:  0.0000
#>   95% CI:   [4.783832, 5.203838]
```

### ANCOVA

Adjusts for covariates via Y ~ 1 + Z + X + Z\*X with HC0 robust SEs.
Achieves efficiency gains when covariates predict the outcome:

``` r

est_via_ancova(exp_obj)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      5.049752
#>   SE:       0.062919
#>   z:        80.2582
#>   p-value:  0.0000
#>   95% CI:   [4.926433, 5.173071]
```

### Stratified Estimator

Weighted average of within-stratum difference-in-means:

``` r

strata <- ifelse(X[,1] > 0, 1, 0)
est_via_strata(exp_obj, strata = strata)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      4.987470
#>   SE:       0.094578
#>   z:        52.7339
#>   p-value:  0.0000
#>   95% CI:   [4.802101, 5.172840]
```

### Fisher Randomization Test

Tests the sharp null of no treatment effect by comparing the observed DM
statistic to its randomization distribution:

``` r

pval <- test_via_fisher(exp_obj, n_draws = 500)
cat("Fisher p-value:", pval, "\n")
#> Fisher p-value: 0
```

## Network Interference

When units within clusters can affect each other’s outcomes, standard
estimators are biased.
[`clustered()`](https://xiangao.github.io/RCausalModel/reference/clustered.md)
implements the framework of Qu et al. (2021).

### Generate clustered data

``` r

cdata <- generate_fixed_cluster(
  clusters = 200,
  group_struct = c(2, 3),
  tau = 1,
  gamma = c(0.5, 0.1)
)
```

### Clustered AIPW

``` r

cl_obj <- clustered(
  Y = cdata$Y,
  Z = cdata$Z,
  X = cdata$X,
  cluster_labels = cdata$cluster_labels,
  group_labels = cdata$group_labels,
  ingroup_labels = cdata$ingroup_labels,
  n_matches = 50
)
result <- est_via_aipw(cl_obj)
#> Warning in predict.lm(model, newdata = df): prediction from rank-deficient fit;
#> attr(*, "non-estim") has doubtful cases
```

The result is a list of group-specific estimates. Each element contains
`beta_g` (treatment effects at each neighborhood configuration) and `se`
(standard errors):

``` r

for (j in seq_along(result)) {
  cat(sprintf("Group %d:\n", j - 1))
  valid <- !is.nan(result[[j]]$beta_g) & result[[j]]$se > 0
  df <- data.frame(
    g_index = which(valid) - 1,
    beta_g = round(result[[j]]$beta_g[valid], 3),
    se = round(result[[j]]$se[valid], 3)
  )
  print(df, row.names = FALSE)
  cat("\n")
}
#> Group 0:
#>  g_index beta_g    se
#>        0  1.261 0.448
#>        1  1.852 0.442
#>        3  1.043 0.197
#>        4  1.811 0.221
#>        6  1.161 0.221
#>        7  1.976 0.261
#>        9 -5.919 0.354
#>       10  0.829 0.409
#> 
#> Group 1:
#>  g_index beta_g    se
#>        0  1.362 0.354
#>        1  1.777 0.247
#>        2  1.783 0.363
#>        3  0.566 0.290
#>        4  1.659 0.194
#>        5  2.355 0.274
#>        6  1.576 0.352
#>        7  1.829 0.225
#>        8  2.715 0.359
```

### Clustered IPW

``` r

result_ipw <- est_via_ipw(cl_obj)
```

``` r

cat("Clustered IPW beta_g (group 0):\n")
#> Clustered IPW beta_g (group 0):
valid <- !is.nan(result_ipw[[1]]$beta_g) & result_ipw[[1]]$se > 0
round(result_ipw[[1]]$beta_g[valid], 3)
#> [1] 1.248 1.743 1.012 1.750 1.348 2.076 1.628 2.124
```

## Custom Models

You can plug in custom ML models by providing a list with `fit()` and
[`predict()`](https://rdrr.io/r/stats/predict.html) (or
`predict_proba()` for classifiers).

### Random forest AIPW

``` r

rf_out <- cm_random_forest_regressor(ntree = 100)
rf_ps <- cm_random_forest_classifier(ntree = 100)
est_via_aipw(obs, outcome_model = rf_out, propensity_model = rf_ps)
#> Warning in fix_propensity(ps, obj$eps): Propensity scores had 115 values of 0
#> or 1, trimmed to [0.0001, 0.9999].
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      10.044085
#>   SE:       0.026667
#>   z:        376.6465
#>   p-value:  0.0000
#>   95% CI:   [9.991819, 10.096352]
```

### Custom model protocol

Any model wrapper that follows this protocol works:

``` r

my_ridge <- list(
  fit = function(X, y) {
    X <- as.matrix(X)
    lambda <- 1.0
    p <- ncol(X)
    beta <- solve(crossprod(X) + lambda * diag(p), crossprod(X, y))
    list(beta = beta)
  },
  predict = function(model, X) {
    as.numeric(as.matrix(X) %*% model$beta)
  }
)

est_via_aipw(obs, outcome_model = my_ridge)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      10.071861
#>   SE:       0.523981
#>   z:        19.2218
#>   p-value:  0.0000
#>   95% CI:   [9.044877, 11.098845]
```
