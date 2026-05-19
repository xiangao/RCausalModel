# IPW Estimator

Inverse Probability Weighting estimator for ATE.

## Usage

``` r
est_via_ipw(obj, ...)

# S3 method for class 'observational'
est_via_ipw(
  obj,
  propensity_model = NULL,
  propensity = NULL,
  normalize = TRUE,
  ...
)
```

## Arguments

- obj:

  An `observational` or `clustered` object.

- ...:

  Additional arguments passed to methods.

- propensity_model:

  A model wrapper with `fit` and `predict_proba`. Default uses logistic
  regression.

- propensity:

  Optional numeric vector of pre-estimated propensity scores.

- normalize:

  Logical. Normalize weights (Hajek estimator). Default TRUE.

## Value

A `cm_result` object.
