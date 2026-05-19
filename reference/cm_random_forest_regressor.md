# Random Forest Regressor Wrapper

Creates a model wrapper around
[`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
for regression. Requires the randomForest package.

## Usage

``` r
cm_random_forest_regressor(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html).

## Value

A list with `fit` and `predict` functions.
