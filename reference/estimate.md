# Default Estimation

Uses AIPW for binary treatment and DML for continuous treatment.

## Usage

``` r
# S3 method for class 'experimental'
estimate(obj, ...)

estimate(obj, ...)

# S3 method for class 'observational'
estimate(obj, ...)
```

## Arguments

- obj:

  An `observational` or `experimental` object.

- ...:

  Additional arguments passed to the selected estimator.

## Value

A `cm_result` object.
