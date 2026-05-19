# Multinomial Logistic Regression Model Wrapper

Creates a model wrapper around
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) for
multiclass classification (e.g., neighborhood treatment propensity).

## Usage

``` r
cm_multi_logistic()
```

## Value

A list with `fit`, `predict_proba`, and `insample_proba` functions.
