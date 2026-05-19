# Matching Estimator

k-nearest-neighbor matching estimator for ATE, with optional bias
adjustment.

## Usage

``` r
est_via_matching(obj, ...)

# S3 method for class 'observational'
est_via_matching(
  obj,
  num_matches = 1,
  num_matches_for_var = NULL,
  bias_adj = FALSE,
  ...
)
```

## Arguments

- obj:

  An `observational` object.

- ...:

  Additional arguments passed to methods.

- num_matches:

  Number of matches per unit (M). Default 1.

- num_matches_for_var:

  Number of matches for variance estimation (J). Default equals
  `num_matches`.

- bias_adj:

  Logical. Apply bias correction using OLS. Default FALSE.

## Value

A `cm_result` object.
