# Compute Rerandomization Balance Criterion

Mahalanobis-type balance criterion from Morgan & Rubin (2015). \$\$m =
\frac{n_t n_c}{n} (\bar{X}\_t - \bar{X}\_c)' \Sigma^{-1} (\bar{X}\_t -
\bar{X}\_c)\$\$

## Usage

``` r
get_balance(Z, X)
```

## Arguments

- Z:

  Numeric vector. Treatment indicator.

- X:

  Numeric matrix. Covariates.

## Value

Numeric scalar. Balance criterion value.
