# Fisher Randomization Test

Computes a p-value by comparing the observed test statistic (DM) to the
distribution under random treatment reassignment.

## Usage

``` r
test_via_fisher(obj, ...)

# S3 method for class 'experimental'
test_via_fisher(obj, n_draws = 1000, ...)
```

## Arguments

- obj:

  An `experimental` object.

- ...:

  Additional arguments passed to methods.

- n_draws:

  Number of random draws. Default 1000.

## Value

Numeric p-value.
