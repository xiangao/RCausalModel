# Create an Observational Study Object

Constructs an observational study data object for estimating causal
effects from non-experimental data.

## Usage

``` r
observational(Y, Z, X)
```

## Arguments

- Y:

  Numeric vector. Outcome variable.

- Z:

  Numeric vector. Treatment variable (binary or continuous).

- X:

  Numeric matrix or vector. Covariates.

## Value

An `observational` object (S3 list).

## Examples

``` r
data <- generate_data(N = 1000, k = 2, tau = 5)
obs <- observational(data$Y, data$Z, data$X)
est_via_aipw(obs)
#> Causal Model Estimation Result
#> ------------------------------
#>   ATE:      4.987031
#>   SE:       0.076564
#>   z:        65.1358
#>   p-value:  0.0000
#>   95% CI:   [4.836969, 5.137093]
```
