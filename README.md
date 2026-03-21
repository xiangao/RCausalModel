# CausalModel

R package for causal inference estimation in observational studies, randomized experiments, and network interference settings. Port of the [Python CausalModel package](https://github.com/freshtaste/CausalModel) (Qu et al. 2021).

## Installation

```r
# install.packages("devtools")
devtools::install_github("username/RCausalModel")
```

## Observational Studies

```r
library(CausalModel)

# Generate synthetic data with known ATE = 10
data <- generate_data(N = 5000, k = 2, tau = 10)
obs <- observational(data$Y, data$Z, data$X)

# Estimators
est_via_ols(obs)        # OLS with HC0 robust SE
est_via_ipw(obs)        # Inverse probability weighting
est_via_aipw(obs)       # Augmented IPW (doubly robust)
est_via_matching(obs)   # Nearest-neighbor matching
estimate(obs)           # Default (AIPW for binary Z)

# Continuous treatment -> Double ML
data_c <- generate_data_continuous(N = 5000, k = 2, tau = 10)
obs_c <- observational(data_c$Y, data_c$Z, data_c$X)
est_via_dml(obs_c)
```

## Randomized Experiments

```r
exp_obj <- experimental(Y, Z, X)

est_via_dm(exp_obj)                    # Difference in means
est_via_strata(exp_obj, strata = s)    # Stratified estimator
est_via_ancova(exp_obj)                # ANCOVA with interactions
test_via_fisher(exp_obj, n_draws = 1000)  # Fisher randomization test
```

## Experimental Designs

```r
design <- crd(treated_ratio = 0.5)        # Completely randomized
design <- bernoulli(treated_prob = 0.5)   # Bernoulli
Z <- draw(design, n = 100)               # Draw treatment vector

# Rerandomization for covariate balance
design <- crd(treated_ratio = 0.5, covariate = X, balance = TRUE, eps = 0.1)
```

## Network Interference

```r
# Generate clustered data with spillover effects
data <- generate_fixed_cluster(clusters = 5000, group_struct = c(2, 3, 4),
                                k = 2, tau = 1, gamma = c(0.5, 0.1, 0.2))

cl <- clustered(data$Y, data$Z, data$X,
                data$cluster_labels, data$group_labels,
                data$ingroup_labels)

# Estimate treatment effects under interference
result <- est_via_ipw(cl)    # IPW under partial interference
result <- est_via_aipw(cl)   # AIPW under partial interference
```

## Custom Models

Supply your own outcome/propensity models:

```r
# Built-in wrappers
est_via_aipw(obs, outcome_model = cm_ols(), propensity_model = cm_logistic())

# Random forest (requires randomForest package)
est_via_aipw(obs,
  outcome_model = cm_random_forest_regressor(),
  propensity_model = cm_random_forest_classifier())
```

## Methods

| Estimator | Function | Setting |
|-----------|----------|---------|
| OLS | `est_via_ols()` | Observational |
| IPW | `est_via_ipw()` | Observational / Interference |
| AIPW | `est_via_aipw()` | Observational / Interference |
| Matching | `est_via_matching()` | Observational |
| Double ML | `est_via_dml()` | Observational (continuous Z) |
| Difference in Means | `est_via_dm()` | Experimental |
| Stratified | `est_via_strata()` | Experimental |
| ANCOVA | `est_via_ancova()` | Experimental |
| Fisher Test | `test_via_fisher()` | Experimental |

## Vignettes

| Vignette | Description |
|----------|-------------|
| `quickstart` | Getting started — basic usage of all estimators |
| `simulations` | Monte Carlo validation across all settings |
| `observational` | Observational estimators with per-estimator simulations (mirrors Python `observational.ipynb`) |
| `experimental` | Design objects and experimental estimators (mirrors Python `experimental.ipynb`) |
| `interference` | Clustered interference with normality and coverage simulations (mirrors Python `Interference.ipynb`) |

## References

Qu, Z., McNamara, R., Jia, R., & Imbens, G. (2021). CausalModel: A Python Package for Causal Inference.
