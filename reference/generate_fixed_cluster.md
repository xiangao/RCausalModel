# Generate Fixed-Size Cluster Data with Interference

Generates clustered data with known direct effect `tau` and spillover
effects `gamma`.

## Usage

``` r
generate_fixed_cluster(
  clusters = 10000,
  group_struct = c(2, 3, 4),
  k = 2,
  tau = 1,
  gamma = NULL,
  label_start = 0L
)
```

## Arguments

- clusters:

  Integer. Number of clusters. Default 10000.

- group_struct:

  Integer vector. Number of units per group within each cluster. E.g.,
  `c(2, 3, 4)` means 2 in group 0, 3 in group 1, 4 in group 2.

- k:

  Integer. Number of individual covariates. Default 2.

- tau:

  Numeric. Direct treatment effect. Default 1.

- gamma:

  Numeric vector or array. Spillover effect parameters. Default is 0.1
  for each group.

- label_start:

  Integer. Starting cluster label. Default 0.

## Value

A list with components: `Y`, `Z`, `X`, `cluster_labels`, `group_labels`,
`ingroup_labels`, `G`, `Xc`.
