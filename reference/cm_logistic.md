# Logistic Regression Model Wrapper

Creates a model wrapper around `glm(family = binomial)` for binary
classification / propensity score estimation.

## Usage

``` r
cm_logistic()
```

## Value

A list with `fit`, `predict_proba`, and `insample_proba` functions.
