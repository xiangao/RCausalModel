# Create a Potential Outcome Data Container

Validates and stores (Y, Z, X) data, splitting into treated/control
groups.

## Usage

``` r
po_data(Y, Z, X = NULL)
```

## Arguments

- Y:

  Numeric vector. Outcome variable.

- Z:

  Numeric vector. Treatment indicator (0/1 for binary, continuous
  allowed).

- X:

  Numeric matrix or vector. Covariates. Vectors are converted to
  single-column matrices.

## Value

A `po_data` object (S3 list).
