# Generate Clustered Data with Varying Cluster Sizes

Combines multiple fixed-cluster datasets with different group
structures.

## Usage

``` r
generate_clustered_data(
  clusters_list = c(5000, 5000, 2000),
  group_struct_list = list(c(2, 3, 4), c(3, 4, 5), c(4, 5, 6)),
  tau = 1,
  gamma = c(5, 0, 1)
)
```

## Arguments

- clusters_list:

  Integer vector. Number of clusters for each structure.

- group_struct_list:

  List of integer vectors. Group structures.

- tau:

  Numeric. Direct treatment effect. Default 1.

- gamma:

  Numeric vector. Spillover parameters.

## Value

A list with concatenated `Y`, `Z`, `X`, `cluster_labels`,
`group_labels`, `ingroup_labels`, `G`, `Xc`.
