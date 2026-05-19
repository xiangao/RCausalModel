# OLS Estimator

Estimates ATE via linear regression of Y on (Z, X) with HC0
heteroskedasticity-robust standard errors.

## Usage

``` r
est_via_ols(obj, ...)

# S3 method for class 'observational'
est_via_ols(obj, ...)
```

## Arguments

- obj:

  An `observational` object.

- ...:

  Additional arguments (currently unused).

## Value

A `cm_result` object.
