# Create an Experimental Study Object

Constructs an experimental study data object for estimating causal
effects from randomized experiments.

## Usage

``` r
experimental(Y, Z, X = NULL, design = NULL)
```

## Arguments

- Y:

  Numeric vector. Outcome variable.

- Z:

  Numeric vector. Binary treatment indicator (0/1).

- X:

  Numeric matrix, vector, or NULL. Covariates (optional).

- design:

  A design object (e.g., from
  [`crd()`](https://xiangao.github.io/RCausalModel/reference/crd.md) or
  [`bernoulli()`](https://xiangao.github.io/RCausalModel/reference/bernoulli.md)).
  Default is CRD.

## Value

An `experimental` object (S3 list).
