#' Generate Observational Data
#'
#' Generates synthetic observational data (Y, Z, X) with known treatment effect
#' \code{tau}. The propensity score is a sigmoid function of the covariates.
#'
#' @param N Integer. Number of observations. Default 10000.
#' @param k Integer. Number of covariates. Default 2.
#' @param tau Numeric. True average treatment effect. Default 10.
#' @return A list with components \code{Y}, \code{Z}, and \code{X}.
#' @export
#' @examples
#' data <- generate_data(N = 1000, k = 2, tau = 5)
#' obs <- observational(data$Y, data$Z, data$X)
#' est_via_aipw(obs)
generate_data <- function(N = 10000, k = 2, tau = 10) {
  X <- matrix(stats::rnorm(N * k), nrow = N, ncol = k)
  beta <- seq(-1, 1, length.out = k)
  prob <- sigmoid(X %*% beta)
  Z <- as.numeric(stats::runif(N) < prob)
  Y <- tau * Z + X %*% beta + stats::rnorm(N)
  list(Y = as.numeric(Y), Z = Z, X = X)
}

#' Generate Observational Data with Continuous Treatment
#'
#' @param N Integer. Number of observations. Default 10000.
#' @param k Integer. Number of covariates. Default 2.
#' @param tau Numeric. True treatment effect. Default 10.
#' @return A list with components \code{Y}, \code{Z}, and \code{X}.
#' @export
generate_data_continuous <- function(N = 10000, k = 2, tau = 10) {
  X <- matrix(stats::rnorm(N * k), nrow = N, ncol = k)
  beta <- seq(-1, 1, length.out = k)
  Z <- as.numeric(X %*% beta) + stats::rnorm(N)
  Y <- tau * Z + as.numeric(X %*% beta) + stats::rnorm(N)
  list(Y = Y, Z = Z, X = X)
}

#' Generate Fixed-Size Cluster Data with Interference
#'
#' Generates clustered data with known direct effect \code{tau} and spillover
#' effects \code{gamma}.
#'
#' @param clusters Integer. Number of clusters. Default 10000.
#' @param group_struct Integer vector. Number of units per group within each
#'   cluster. E.g., \code{c(2, 3, 4)} means 2 in group 0, 3 in group 1, 4 in
#'   group 2.
#' @param k Integer. Number of individual covariates. Default 2.
#' @param tau Numeric. Direct treatment effect. Default 1.
#' @param gamma Numeric vector or array. Spillover effect parameters. Default
#'   is 0.1 for each group.
#' @param label_start Integer. Starting cluster label. Default 0.
#' @return A list with components: \code{Y}, \code{Z}, \code{X},
#'   \code{cluster_labels}, \code{group_labels}, \code{ingroup_labels},
#'   \code{G}, \code{Xc}.
#' @export
generate_fixed_cluster <- function(clusters = 10000, group_struct = c(2, 3, 4),
                                    k = 2, tau = 1.0, gamma = NULL,
                                    label_start = 0L) {
  group_struct <- as.integer(group_struct)
  ngroup <- length(group_struct)

  # Build treatment effect array beta(g) on grid
  grid_dims <- group_struct + 1L
  G_count <- prod(grid_dims)

  # Create meshgrid indices
  grid_list <- lapply(group_struct, function(gs) 0:gs)
  grid <- as.matrix(expand.grid(grid_list))  # G_count x ngroup

  if (is.null(gamma)) {
    gamma_marginal <- rep(0.1, ngroup)
    beta <- tau + as.numeric(grid %*% gamma_marginal)
  } else if (is.vector(gamma) && length(gamma) == ngroup) {
    beta <- tau + as.numeric(grid %*% gamma)
  } else if (is.array(gamma) || is.matrix(gamma)) {
    beta <- tau + as.numeric(gamma)
  } else {
    stop("gamma shape mismatch")
  }

  # Beta as array indexed by G
  beta_arr <- array(beta, dim = grid_dims)

  nunit_per_cluster <- sum(group_struct)
  units <- clusters * nunit_per_cluster

  # Labels
  cluster_labels <- label_start + rep(0:(clusters - 1), each = nunit_per_cluster)
  group_labels_wc <- rep(0:(ngroup - 1), times = group_struct)
  group_labels <- rep(group_labels_wc, clusters)
  ingroup_labels_wc <- unlist(lapply(group_struct, function(g) 0:(g - 1)))
  ingroup_labels <- rep(ingroup_labels_wc, clusters)

  # Covariates: X ~ N(0, 0.01*I) per unit
  X <- 0.1 * matrix(stats::rnorm(units * k), nrow = units, ncol = k)
  X_3d <- array(X, dim = c(nunit_per_cluster, clusters, k))
  X_3d <- aperm(X_3d, c(2, 1, 3))  # clusters x nunit x k

  # Cluster mean
  Xcmean <- apply(X_3d, c(1, 3), mean)  # clusters x k

  # Augmented covariates: [X, Xcmean]
  Xc <- matrix(0, nrow = units, ncol = 2 * k)
  Xc[, 1:k] <- X
  Xcmean_rep <- Xcmean[rep(1:clusters, each = nunit_per_cluster), ]
  Xc[, (k + 1):(2 * k)] <- Xcmean_rep

  # Treatment assignment
  beta_xc <- seq(-1, 1, length.out = 2 * k)
  prop_idv <- sigmoid(Xc %*% beta_xc)
  Z <- as.integer(stats::runif(units) < prop_idv)

  # Compute G (number of treated neighbours, excluding self)
  Z_onehot <- array(0L, dim = c(clusters, nunit_per_cluster, ngroup))
  cl_offset <- cluster_labels - label_start
  unit_in_cl <- rep(1:nunit_per_cluster, clusters)
  for (i in seq_len(units)) {
    Z_onehot[cl_offset[i] + 1, unit_in_cl[i], group_labels[i] + 1] <- Z[i]
  }

  G_plus <- array(0L, dim = c(clusters, nunit_per_cluster, ngroup))
  for (g in seq_len(ngroup)) {
    total_g <- rowSums(matrix(Z_onehot[, , g], nrow = clusters))
    G_plus[, , g] <- matrix(total_g, nrow = clusters, ncol = nunit_per_cluster)
  }
  G_stack <- G_plus - Z_onehot
  G <- matrix(aperm(G_stack, c(2, 1, 3)), nrow = units, ncol = ngroup)

  # Outcome
  epsilon <- stats::rnorm(units)
  # beta(G) for each unit: index into beta_arr
  te <- numeric(units)
  for (i in seq_len(units)) {
    idx <- as.list(G[i, ] + 1L)  # 1-indexed for R arrays
    te[i] <- Z[i] * do.call("[", c(list(beta_arr), idx))
  }

  Y <- (1 + Z) * as.numeric(Xc %*% beta_xc) + te + epsilon

  # Shuffle
  perm <- sample(units)
  list(
    Y = Y[perm],
    Z = Z[perm],
    X = X[perm, ],
    cluster_labels = cluster_labels[perm],
    group_labels = group_labels[perm],
    ingroup_labels = ingroup_labels[perm],
    G = G[perm, ],
    Xc = Xc[perm, ]
  )
}

#' Generate Clustered Data with Varying Cluster Sizes
#'
#' Combines multiple fixed-cluster datasets with different group structures.
#'
#' @param clusters_list Integer vector. Number of clusters for each structure.
#' @param group_struct_list List of integer vectors. Group structures.
#' @param tau Numeric. Direct treatment effect. Default 1.
#' @param gamma Numeric vector. Spillover parameters.
#' @return A list with concatenated \code{Y}, \code{Z}, \code{X},
#'   \code{cluster_labels}, \code{group_labels}, \code{ingroup_labels},
#'   \code{G}, \code{Xc}.
#' @export
generate_clustered_data <- function(clusters_list = c(5000, 5000, 2000),
                                     group_struct_list = list(c(2, 3, 4),
                                                               c(3, 4, 5),
                                                               c(4, 5, 6)),
                                     tau = 1.0,
                                     gamma = c(5, 0, 1)) {
  label_starts <- c(0, cumsum(clusters_list[-length(clusters_list)]))

  datasets <- mapply(function(cl, gs, ls) {
    generate_fixed_cluster(clusters = cl, group_struct = gs,
                           tau = tau, gamma = gamma, label_start = ls)
  }, clusters_list, group_struct_list, label_starts, SIMPLIFY = FALSE)

  # Concatenate
  list(
    Y = unlist(lapply(datasets, `[[`, "Y")),
    Z = unlist(lapply(datasets, `[[`, "Z")),
    X = do.call(rbind, lapply(datasets, `[[`, "X")),
    cluster_labels = unlist(lapply(datasets, `[[`, "cluster_labels")),
    group_labels = unlist(lapply(datasets, `[[`, "group_labels")),
    ingroup_labels = unlist(lapply(datasets, `[[`, "ingroup_labels")),
    G = do.call(rbind, lapply(datasets, `[[`, "G")),
    Xc = do.call(rbind, lapply(datasets, `[[`, "Xc"))
  )
}


#' @keywords internal
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}
