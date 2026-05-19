# Completely Randomized Design

Creates a CRD design object for random treatment assignment. Optionally
supports rerandomization for covariate balance.

## Usage

``` r
crd(
  treated_ratio = 0.5,
  covariate = NULL,
  balance = FALSE,
  eps = 0.1,
  max_iter = 1000
)
```

## Arguments

- treated_ratio:

  Numeric. Fraction of units assigned to treatment. Default 0.5.

- covariate:

  Numeric matrix. Covariates for balance checking (optional).

- balance:

  Logical. Use rerandomization for covariate balance. Default FALSE.

- eps:

  Numeric. Balance criterion threshold. Default 0.1.

- max_iter:

  Integer. Maximum rerandomization attempts. Default 1000.

## Value

A `crd` design object.
