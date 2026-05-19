# Create a Causal Model Result

Stores the output of a causal inference estimator: point estimate,
standard error, z-statistic, p-value, and 95% confidence interval.

## Usage

``` r
cm_result(ate, se, alpha = 0.05)
```

## Arguments

- ate:

  Numeric. Average treatment effect estimate.

- se:

  Numeric. Standard error of the ATE.

- alpha:

  Numeric. Significance level for CI (default 0.05).

## Value

A `cm_result` object (S3 list).
