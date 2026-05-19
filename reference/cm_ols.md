# OLS Model Wrapper

Creates a model wrapper around [`lm()`](https://rdrr.io/r/stats/lm.html)
that conforms to the CausalModel model protocol: a list with `fit(X, y)`
and `predict(model, X)` functions.

## Usage

``` r
cm_ols()
```

## Value

A list with `fit` and `predict` functions.
