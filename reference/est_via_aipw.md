# AIPW (Doubly Robust) Estimator

Augmented IPW estimator. Consistent if either the outcome model or the
propensity model is correctly specified.

## Usage

``` r
est_via_aipw(obj, ...)

# S3 method for class 'observational'
est_via_aipw(
  obj,
  outcome_model = NULL,
  propensity_model = NULL,
  treated_pred = NULL,
  control_pred = NULL,
  propensity = NULL,
  ...
)
```

## Arguments

- obj:

  An `observational` or `clustered` object.

- ...:

  Additional arguments passed to methods.

- outcome_model:

  A model wrapper with `fit` and `predict`. Default uses OLS.

- propensity_model:

  A model wrapper with `fit` and `predict_proba`. Default uses logistic
  regression.

- treated_pred:

  Optional pre-computed treated outcome predictions.

- control_pred:

  Optional pre-computed control outcome predictions.

- propensity:

  Optional pre-computed propensity scores.

## Value

A `cm_result` object.
