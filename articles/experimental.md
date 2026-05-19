# Experimental Designs and Estimators

## Overview

This vignette mirrors the Python `experimental.ipynb` test notebook,
demonstrating the design objects and experimental estimators in
CausalModel.

## Design Objects

CausalModel provides two randomization designs: completely randomized
design (CRD) and Bernoulli randomization.

### Completely Randomized Design

``` r

d <- crd()
draw(d, 10)
#>  [1] 1 1 0 1 1 0 0 0 0 1
```

### CRD with Covariate Balance

Rerandomization to ensure covariate balance (Morgan & Rubin, 2012). The
design draws assignments until Mahalanobis balance falls below a
threshold:

``` r

X <- matrix(rnorm(100 * 2), ncol = 2)
d_bal <- crd(covariate = X, balance = TRUE)
Z_bal <- draw(d_bal, 100)
cat("Treatment assigned:", sum(Z_bal), "of", length(Z_bal), "\n")
#> Treatment assigned: 50 of 100
```

### Bernoulli Randomization

Each unit is independently assigned to treatment with a fixed
probability:

``` r

d_bnli <- bernoulli()
draw(d_bnli, 10)
#>  [1] 0 0 1 0 0 0 1 0 0 1
```

## Experimental Estimators

### Setup

Generate experimental data with a known treatment effect $`\tau = 10`$:

``` r

data <- generate_data(N = 10000, k = 2, tau = 10)
exp_obj <- experimental(data$Y, data$Z, data$X)
```

### Difference-in-Means

The default estimator under a CRD:

``` r

est_via_dm(exp_obj)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      11.483717
#>   SE:       0.031312
#>   z:        366.7547
#>   p-value:  0.0000
#>   95% CI:   [11.422347, 11.545087]
```

### Stratified Estimator

Weighted average of within-stratum difference-in-means. Strata are
defined by the user — here we use a simple median split on the first
covariate:

``` r

strata <- as.integer(data$X[, 1] > median(data$X[, 1]))
est_via_strata(exp_obj, strata = strata)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      11.088742
#>   SE:       0.030043
#>   z:        369.0955
#>   p-value:  0.0000
#>   95% CI:   [11.029859, 11.147626]
```

With random strata (less informative):

``` r

strata_rand <- sample(0:2, length(data$Y), replace = TRUE)
est_via_strata(exp_obj, strata = strata_rand)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      11.484102
#>   SE:       0.031298
#>   z:        366.9291
#>   p-value:  0.0000
#>   95% CI:   [11.422759, 11.545444]
```

### ANCOVA

Adjusts for covariates via $`Y \sim 1 + Z + X + Z \cdot X`$ with HC0
robust standard errors. Achieves efficiency gains when covariates are
predictive:

``` r

est_via_ancova(exp_obj)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      10.029845
#>   SE:       0.023435
#>   z:        427.9815
#>   p-value:  0.0000
#>   95% CI:   [9.983913, 10.075777]
```

### Fisher Randomization Test

Tests the sharp null hypothesis of zero individual treatment effects by
comparing the observed test statistic to its randomization distribution:

``` r

pval <- test_via_fisher(exp_obj, n_draws = 1000)
cat("Fisher p-value:", pval, "\n")
#> Fisher p-value: 0
```

With $`\tau = 10`$, the Fisher test rejects decisively (p-value
$`\approx 0`$).

## Comparing Estimators

Under a CRD with predictive covariates, ANCOVA should achieve
substantially smaller standard errors than DM while both remain
unbiased:

``` r

results <- data.frame(
  Estimator = c("DM", "Strata (median split)", "ANCOVA"),
  ATE = NA, SE = NA
)

r <- est_via_dm(exp_obj)
results[1, c("ATE", "SE")] <- c(r$ate, r$se)

r <- est_via_strata(exp_obj, strata = strata)
results[2, c("ATE", "SE")] <- c(r$ate, r$se)

r <- est_via_ancova(exp_obj)
results[3, c("ATE", "SE")] <- c(r$ate, r$se)

knitr::kable(results, digits = 4,
             caption = "Estimator comparison (N=10000, tau=10)")
```

| Estimator             |     ATE |     SE |
|:----------------------|--------:|-------:|
| DM                    | 11.4837 | 0.0313 |
| Strata (median split) | 11.0887 | 0.0300 |
| ANCOVA                | 10.0298 | 0.0234 |

Estimator comparison (N=10000, tau=10) {.table}

## Summary

- **Design objects**
  ([`crd()`](https://xiangao.github.io/RCausalModel/reference/crd.md),
  [`bernoulli()`](https://xiangao.github.io/RCausalModel/reference/bernoulli.md))
  generate treatment assignments via
  [`draw()`](https://xiangao.github.io/RCausalModel/reference/draw.md).
  CRD with `balance = TRUE` ensures covariate balance through
  rerandomization.
- **DM** is the natural estimator under randomization but ignores
  covariate information.
- **Stratified estimation** improves precision when strata are
  informative.
- **ANCOVA** provides the largest efficiency gain when covariates
  predict the outcome.
- **Fisher test** provides exact p-values under the sharp null, without
  asymptotic approximations.
