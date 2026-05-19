# Generate Observational Data with Continuous Treatment

Generate Observational Data with Continuous Treatment

## Usage

``` r
generate_data_continuous(N = 10000, k = 2, tau = 10)
```

## Arguments

- N:

  Integer. Number of observations. Default 10000.

- k:

  Integer. Number of covariates. Default 2.

- tau:

  Numeric. True treatment effect. Default 10.

## Value

A list with components `Y`, `Z`, and `X`.
