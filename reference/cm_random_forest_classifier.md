# Random Forest Classifier Wrapper

Creates a model wrapper around
[`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
for binary classification / propensity score estimation. Requires the
randomForest package.

## Usage

``` r
cm_random_forest_classifier(...)
```

## Arguments

- ...:

  Additional arguments passed to
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html).

## Value

A list with `fit`, `predict_proba`, and `insample_proba` functions.
