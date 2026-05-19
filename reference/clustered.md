# Create a Clustered Interference Study Object

Estimates causal effects under partial/clustered interference, where
units within clusters may affect each other's outcomes.

## Usage

``` r
clustered(
  Y,
  Z,
  X,
  cluster_labels,
  group_labels,
  ingroup_labels,
  cluster_feature = NULL,
  n_moments = 1L,
  prop_idv_model = NULL,
  prop_neigh_model = NULL,
  n_matches = 100L,
  subsampling_match = 2000L,
  categorical_Z = TRUE
)
```

## Arguments

- Y:

  Numeric vector. Outcome variable.

- Z:

  Numeric vector. Treatment indicator (0/1).

- X:

  Numeric matrix. Individual-level covariates.

- cluster_labels:

  Integer vector. Cluster identifier for each unit.

- group_labels:

  Integer vector. Group index within cluster (0-indexed).

- ingroup_labels:

  Integer vector. Position within group (0-indexed).

- cluster_feature:

  Optional numeric matrix. Cluster-level covariates.

- n_moments:

  Integer. Number of moments for covariate aggregation. Default 1.

- prop_idv_model:

  Model wrapper for individual propensity P(Z=1\|X).

- prop_neigh_model:

  Model wrapper for neighborhood propensity P(G\|X).

- n_matches:

  Integer. Number of matches for variance estimation. Default 100.

- subsampling_match:

  Integer. Maximum subsample size for matching. Default 2000.

- categorical_Z:

  Logical. Treat Z as categorical. Default TRUE.

## Value

A `clustered` object.
