# Generate Observational Data

Generates synthetic observational data (Y, Z, X) with known treatment
effect `tau`. The propensity score is a sigmoid function of the
covariates.

## Usage

``` r
generate_data(N = 10000, k = 2, tau = 10)
```

## Arguments

- N:

  Integer. Number of observations. Default 10000.

- k:

  Integer. Number of covariates. Default 2.

- tau:

  Numeric. True average treatment effect. Default 10.

## Value

A list with components `Y`, `Z`, and `X`.

## Examples

``` r
data <- generate_data(N = 1000, k = 2, tau = 5)
obs <- observational(data$Y, data$Z, data$X)
est_via_aipw(obs)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      4.929815
#>   SE:       0.073875
#>   z:        66.7319
#>   p-value:  0.0000
#>   95% CI:   [4.785023, 5.074608]
```
