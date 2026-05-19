# Package index

## Study designs

Constructors for observational, experimental, and clustered-interference
designs.

- [`observational()`](https://xiangao.github.io/RCausalModel/reference/observational.md)
  : Create an Observational Study Object
- [`experimental()`](https://xiangao.github.io/RCausalModel/reference/experimental.md)
  : Create an Experimental Study Object
- [`clustered()`](https://xiangao.github.io/RCausalModel/reference/clustered.md)
  : Create a Clustered Interference Study Object
- [`crd()`](https://xiangao.github.io/RCausalModel/reference/crd.md) :
  Completely Randomized Design
- [`bernoulli()`](https://xiangao.github.io/RCausalModel/reference/bernoulli.md)
  : Bernoulli Randomized Design
- [`draw()`](https://xiangao.github.io/RCausalModel/reference/draw.md) :
  Draw Treatment Assignment
- [`get_params_via_obs()`](https://xiangao.github.io/RCausalModel/reference/get_params_via_obs.md)
  : Infer Design Parameters from Observed Treatment

## Estimators

Estimation methods for causal effects and randomization tests.

- [`estimate()`](https://xiangao.github.io/RCausalModel/reference/estimate.md)
  : Default Estimation
- [`est_via_ols()`](https://xiangao.github.io/RCausalModel/reference/est_via_ols.md)
  : OLS Estimator
- [`est_via_ipw()`](https://xiangao.github.io/RCausalModel/reference/est_via_ipw.md)
  : IPW Estimator
- [`est_via_aipw()`](https://xiangao.github.io/RCausalModel/reference/est_via_aipw.md)
  : AIPW (Doubly Robust) Estimator
- [`est_via_matching()`](https://xiangao.github.io/RCausalModel/reference/est_via_matching.md)
  : Matching Estimator
- [`est_via_dml()`](https://xiangao.github.io/RCausalModel/reference/est_via_dml.md)
  : Double/Debiased Machine Learning Estimator
- [`est_via_dm()`](https://xiangao.github.io/RCausalModel/reference/est_via_dm.md)
  : Difference-in-Means Estimator
- [`est_via_strata()`](https://xiangao.github.io/RCausalModel/reference/est_via_strata.md)
  : Stratified Estimator
- [`est_via_ancova()`](https://xiangao.github.io/RCausalModel/reference/est_via_ancova.md)
  : ANCOVA Estimator
- [`test_via_fisher()`](https://xiangao.github.io/RCausalModel/reference/test_via_fisher.md)
  : Fisher Randomization Test

## Data generation and models

- [`generate_data()`](https://xiangao.github.io/RCausalModel/reference/generate_data.md)
  : Generate Observational Data
- [`generate_data_continuous()`](https://xiangao.github.io/RCausalModel/reference/generate_data_continuous.md)
  : Generate Observational Data with Continuous Treatment
- [`generate_clustered_data()`](https://xiangao.github.io/RCausalModel/reference/generate_clustered_data.md)
  : Generate Clustered Data with Varying Cluster Sizes
- [`generate_fixed_cluster()`](https://xiangao.github.io/RCausalModel/reference/generate_fixed_cluster.md)
  : Generate Fixed-Size Cluster Data with Interference
- [`po_data()`](https://xiangao.github.io/RCausalModel/reference/po_data.md)
  : Create a Potential Outcome Data Container
- [`get_balance()`](https://xiangao.github.io/RCausalModel/reference/get_balance.md)
  : Compute Rerandomization Balance Criterion
- [`cm_ols()`](https://xiangao.github.io/RCausalModel/reference/cm_ols.md)
  : OLS Model Wrapper
- [`cm_logistic()`](https://xiangao.github.io/RCausalModel/reference/cm_logistic.md)
  : Logistic Regression Model Wrapper
- [`cm_multi_logistic()`](https://xiangao.github.io/RCausalModel/reference/cm_multi_logistic.md)
  : Multinomial Logistic Regression Model Wrapper
- [`cm_random_forest_regressor()`](https://xiangao.github.io/RCausalModel/reference/cm_random_forest_regressor.md)
  : Random Forest Regressor Wrapper
- [`cm_random_forest_classifier()`](https://xiangao.github.io/RCausalModel/reference/cm_random_forest_classifier.md)
  : Random Forest Classifier Wrapper
- [`cm_result()`](https://xiangao.github.io/RCausalModel/reference/cm_result.md)
  : Create a Causal Model Result
