# Stratified Estimator

Estimates ATE as a weighted average of within-stratum
difference-in-means.

## Usage

``` r
est_via_strata(obj, ...)

# S3 method for class 'experimental'
est_via_strata(obj, strata, ...)
```

## Arguments

- obj:

  An `experimental` object.

- ...:

  Additional arguments passed to methods.

- strata:

  Integer or factor vector of stratum labels (same length as Y).

## Value

A `cm_result` object.
