# Bernoulli Randomized Design

Creates a Bernoulli design where each unit is independently assigned to
treatment with probability `treated_prob`.

## Usage

``` r
bernoulli(
  treated_prob = 0.5,
  covariate = NULL,
  balance = FALSE,
  eps = 0.1,
  max_iter = 1000
)
```

## Arguments

- treated_prob:

  Numeric. Probability of treatment. Default 0.5.

- covariate:

  Numeric matrix. Covariates for balance checking (optional).

- balance:

  Logical. Use rerandomization. Default FALSE.

- eps:

  Numeric. Balance criterion threshold. Default 0.1.

- max_iter:

  Integer. Maximum rerandomization attempts. Default 1000.

## Value

A `bernoulli` design object.
