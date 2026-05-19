# ANCOVA Estimator

Fisher's ANCOVA with treatment-covariate interactions and HC0 robust
standard errors. Regresses Y ~ 1 + Z + X + Z\*X.

## Usage

``` r
est_via_ancova(obj, ...)

# S3 method for class 'experimental'
est_via_ancova(obj, ...)
```

## Arguments

- obj:

  An `experimental` object.

- ...:

  Additional arguments (currently unused).

## Value

A `cm_result` object.
