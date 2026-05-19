# Double/Debiased Machine Learning Estimator

DML estimator for continuous treatments using cross-fitting.

## Usage

``` r
est_via_dml(obj, ...)

# S3 method for class 'observational'
est_via_dml(
  obj,
  outcome_model = NULL,
  treatment_model = NULL,
  k_folds = 2,
  ...
)
```

## Arguments

- obj:

  An `observational` object.

- ...:

  Additional arguments passed to methods.

- outcome_model:

  Model wrapper for Y ~ X. Default uses OLS.

- treatment_model:

  Model wrapper for Z ~ X. Default uses OLS.

- k_folds:

  Number of cross-fitting folds. Default 2.

## Value

A `cm_result` object.
