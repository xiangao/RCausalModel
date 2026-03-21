#' Create a Clustered Interference Study Object
#'
#' Estimates causal effects under partial/clustered interference, where units
#' within clusters may affect each other's outcomes.
#'
#' @param Y Numeric vector. Outcome variable.
#' @param Z Numeric vector. Treatment indicator (0/1).
#' @param X Numeric matrix. Individual-level covariates.
#' @param cluster_labels Integer vector. Cluster identifier for each unit.
#' @param group_labels Integer vector. Group index within cluster (0-indexed).
#' @param ingroup_labels Integer vector. Position within group (0-indexed).
#' @param cluster_feature Optional numeric matrix. Cluster-level covariates.
#' @param n_moments Integer. Number of moments for covariate aggregation.
#'   Default 1.
#' @param prop_idv_model Model wrapper for individual propensity P(Z=1|X).
#' @param prop_neigh_model Model wrapper for neighborhood propensity P(G|X).
#' @param n_matches Integer. Number of matches for variance estimation.
#'   Default 100.
#' @param subsampling_match Integer. Maximum subsample size for matching.
#'   Default 2000.
#' @param categorical_Z Logical. Treat Z as categorical. Default TRUE.
#' @return A \code{clustered} object.
#' @export
clustered <- function(Y, Z, X, cluster_labels, group_labels, ingroup_labels,
                      cluster_feature = NULL, n_moments = 1L,
                      prop_idv_model = NULL, prop_neigh_model = NULL,
                      n_matches = 100L, subsampling_match = 2000L,
                      categorical_Z = TRUE) {
  if (is.null(prop_idv_model)) prop_idv_model <- cm_logistic()
  if (is.null(prop_neigh_model)) prop_neigh_model <- cm_multi_logistic()

  cluster_labels <- as.integer(cluster_labels)
  group_labels <- as.integer(group_labels)
  ingroup_labels <- as.integer(ingroup_labels)

  cdata <- cluster_data(Y, Z, X, cluster_labels, group_labels, ingroup_labels,
                        cluster_feature, n_moments, categorical_Z)

  obj <- list(
    data = cdata,
    prop_idv_model = prop_idv_model,
    prop_neigh_model = prop_neigh_model,
    n_matches = n_matches,
    subsampling_match = subsampling_match,
    eps = 1e-4
  )
  class(obj) <- c("clustered", "observational")
  obj
}

#' @export
est_via_ipw.clustered <- function(obj, ...) {
  clustered_est(obj, method = "ipw")
}

#' @export
est_via_aipw.clustered <- function(obj, ...) {
  clustered_est(obj, method = "aipw")
}


# ---- Core estimation engine ----

#' @keywords internal
clustered_est <- function(obj, method = "ipw") {
  cdata <- obj$data
  group_structs <- sort(names(cdata$data_by_group_struct))
  n_groups <- length(cdata$first_group_struct)

  total_result <- lapply(seq_len(n_groups), function(j) list())

  for (gs_name in group_structs) {
    group_struct <- cdata$data_by_group_struct[[gs_name]]$group_struct
    sub <- cdata$data_by_group_struct[[gs_name]]
    individuals_in_cluster <- sum(group_struct)
    Mn <- length(sub$Y) / individuals_in_cluster
    pn <- Mn / cdata$M

    result <- est_subsample(obj, sub, group_struct, method)

    for (j in seq_len(n_groups)) {
      taug <- result[[j]]$beta_g
      Vg <- result[[j]]$se^2 * Mn
      total_result[[j]][[gs_name]] <- list(pn = pn, taug = taug, Vg = Vg,
                                           group_struct = group_struct)
    }
  }

  # Aggregate across group structures
  all_gs <- lapply(names(cdata$data_by_group_struct), function(nm) {
    cdata$data_by_group_struct[[nm]]$group_struct
  })
  max_group_struct <- do.call(pmax, all_gs)

  G_count <- prod(max_group_struct + 1)

  ret <- list()
  for (j in seq_len(n_groups)) {
    beta_g <- rep(0, G_count)
    se_g <- rep(0, G_count)

    for (g_enc in seq_len(G_count) - 1) {  # 0-indexed
      g <- decode_G(g_enc, max_group_struct)

      key_vals <- list()
      for (entry in total_result[[j]]) {
        gg <- entry$group_struct
        if (all(gg >= g) && !all(gg == g)) {
          encoded <- encode_G_single(g, gg) + 1  # R 1-indexed
          kv <- c(entry$pn, entry$taug[encoded], entry$Vg[encoded])
          key_vals <- c(key_vals, list(kv))
        }
      }

      if (length(key_vals) == 0) next

      key_vals <- do.call(rbind, key_vals)
      w <- key_vals[, 1] / sum(key_vals[, 1])
      taug_n <- key_vals[, 2]

      invalid <- is.nan(taug_n)
      taug_n[invalid] <- 0
      w[invalid] <- 0
      if (sum(w) > 0) w <- w / sum(w)

      beta_g[g_enc + 1] <- sum(w * taug_n)

      Vg_n <- key_vals[, 3]
      Vg_n[invalid] <- 0
      Vg1 <- sum(Vg_n * w^2 / key_vals[, 1])
      se_g[g_enc + 1] <- sqrt(Vg1 / cdata$M)
    }

    ret[[j]] <- list(beta_g = beta_g, se = se_g,
                     max_group_struct = max_group_struct)
  }

  ret
}


#' @keywords internal
est_subsample <- function(obj, sub, group_struct, method = "ipw") {
  Y <- sub$Y
  Z <- sub$Z
  G <- sub$G
  Xc <- sub$X
  cluster_labels <- sub$cluster_labels
  group_labels <- sub$group_labels
  ingroup_labels <- sub$ingroup_labels

  M <- length(unique(cluster_labels))
  N <- length(Y)

  G_encoded <- encode_G_mat(G, group_struct)

  # Standardize covariates
  sd_Xc <- apply(Xc, 2, function(x) sqrt(mean((x - mean(x))^2)))
  sd_Xc[sd_Xc == 0] <- 1
  Xc_s <- sweep(Xc, 2, sd_Xc, "/")
  X_g <- cbind(Xc_s, group_labels)

  # Estimate propensities
  fitted_idv <- obj$prop_idv_model$fit(X_g, Z)
  prop_idv <- obj$prop_idv_model$insample_proba(fitted_idv)

  fitted_neigh <- obj$prop_neigh_model$fit(X_g, G_encoded)
  prop_neigh <- obj$prop_neigh_model$insample_proba(fitted_neigh)
  if (is.null(dim(prop_neigh))) {
    prop_neigh <- cbind(1 - prop_neigh, prop_neigh)
  }

  G_count <- prod(group_struct + 1)
  n_groups <- length(group_struct)

  result <- lapply(seq_len(n_groups), function(j) {
    list(beta_g = rep(0, G_count), se = rep(0, G_count))
  })

  ols_model <- cm_ols()

  for (g_enc in seq_len(G_count) - 1) {  # 0-indexed
    g <- decode_G(g_enc, group_struct)

    # Check if we have treated and control with G == g
    g_match <- apply(G, 1, function(row) all(row == g))
    mask1 <- g_match & (Z == 1)
    mask0 <- g_match & (Z == 0)

    if (!any(mask1) || !any(mask0)) {
      for (j in seq_len(n_groups)) {
        result[[j]]$beta_g[g_enc + 1] <- NaN
      }
      # Variance
      for (j in seq_len(n_groups)) {
        result[[j]]$se[g_enc + 1] <- NaN
      }
      next
    }

    if (method == "aipw") {
      fitted_mu1 <- ols_model$fit(Xc[mask1, , drop = FALSE], Y[mask1])
      mu1g <- ols_model$predict(fitted_mu1, Xc)
      fitted_mu0 <- ols_model$fit(Xc[mask0, , drop = FALSE], Y[mask0])
      mu0g <- ols_model$predict(fitted_mu0, Xc)
    }

    for (j in seq_len(n_groups)) {
      mask_j <- group_labels == (j - 1)  # 0-indexed group labels
      gg <- g
      gg[j] <- gg[j] + 1

      if (any(gg > group_struct)) {
        result[[j]]$beta_g[g_enc + 1] <- NaN
      } else if (method == "ipw") {
        result[[j]]$beta_g[g_enc + 1] <- ipw_formula(
          Y[mask_j], Z[mask_j], G[mask_j, , drop = FALSE],
          prop_idv[mask_j], prop_neigh[mask_j, , drop = FALSE],
          g, g_enc)
      } else if (method == "aipw") {
        result[[j]]$beta_g[g_enc + 1] <- aipw_formula(
          Y[mask_j], Z[mask_j], G[mask_j, , drop = FALSE],
          prop_idv[mask_j], prop_neigh[mask_j, , drop = FALSE],
          g, g_enc, mu1g[mask_j], mu0g[mask_j])
      }
    }

    # Variance estimation for each group
    for (j in seq_len(n_groups)) {
      sub_idx <- which(group_labels == (j - 1))
      if (obj$subsampling_match < length(sub_idx)) {
        sub_idx <- sample(sub_idx, obj$subsampling_match, replace = FALSE)
      }

      idx_g <- apply(G[sub_idx, , drop = FALSE], 1, function(row) all(row == g))
      if (!any(idx_g)) {
        result[[j]]$se[g_enc + 1] <- NaN
        next
      }

      Vg <- variance_via_matching(
        obj, Y[sub_idx], Z[sub_idx], Xc[sub_idx, , drop = FALSE],
        ingroup_labels[sub_idx],
        prop_idv[sub_idx], prop_neigh[sub_idx, g_enc + 1],
        group_struct[j], idx_g)
      result[[j]]$se[g_enc + 1] <- sqrt(Vg / M)
    }
  }

  result
}


# ---- IPW and AIPW formulas ----

#' @keywords internal
ipw_formula <- function(Y, Z, G, prop_idv, prop_neigh, g, g_enc) {
  N <- length(Y)
  g_match <- apply(G, 1, function(row) all(row == g))
  w1 <- g_match * Z / (prop_neigh[, g_enc + 1] * prop_idv)
  w0 <- g_match * (1 - Z) / (prop_neigh[, g_enc + 1] * (1 - prop_idv))
  arr <- Y * w1 / (sum(w1) / N) - Y * w0 / (sum(w0) / N)
  mean(arr)
}

#' @keywords internal
aipw_formula <- function(Y, Z, G, prop_idv, prop_neigh, g, g_enc,
                         mu1g, mu0g) {
  N <- length(Y)
  g_match <- apply(G, 1, function(row) all(row == g))
  w1 <- g_match * Z / (prop_neigh[, g_enc + 1] * prop_idv)
  w0 <- g_match * (1 - Z) / (prop_neigh[, g_enc + 1] * (1 - prop_idv))
  arr <- (Y - mu1g) * w1 / (sum(w1) / N) -
    (Y - mu0g) * w0 / (sum(w0) / N) +
    mu1g - mu0g
  mean(arr)
}


# ---- Variance via matching ----

#' @keywords internal
variance_via_matching <- function(obj, Y, Z, Xc, ingroup_labels, q, pg,
                                   size, idx_g) {
  idx1 <- Z == 1
  Y1_g <- Y[idx1 & idx_g]
  Y0_g <- Y[(!idx1) & idx_g]
  Xc1_g <- Xc[idx1 & idx_g, , drop = FALSE]
  Xc0_g <- Xc[(!idx1) & idx_g, , drop = FALSE]

  if (nrow(Xc1_g) < 2 || nrow(Xc0_g) < 2) {
    return(NaN)
  }

  # Standardize
  sd_all <- apply(Xc, 2, function(x) sqrt(mean((x - mean(x))^2)))
  sd_all[sd_all == 0] <- 1
  Xc_s <- sweep(Xc, 2, sd_all, "/")

  sd_1g <- apply(Xc1_g, 2, function(x) sqrt(mean((x - mean(x))^2)))
  sd_1g[sd_1g == 0] <- 1
  X1 <- sweep(Xc1_g, 2, sd_1g, "/")

  sd_0g <- apply(Xc0_g, 2, function(x) sqrt(mean((x - mean(x))^2)))
  sd_0g[sd_0g == 0] <- 1
  X0 <- sweep(Xc0_g, 2, sd_0g, "/")

  # Match all units to treated-G and control-G units
  match_t <- mat_match_mat(Xc_s, X1, obj$n_matches)
  match_c <- mat_match_mat(Xc_s, X0, obj$n_matches)

  # Variance of matched outcomes
  var_t <- apply(matrix(Y1_g[match_t], nrow = nrow(Xc_s)), 1, stats::var)
  var_c <- apply(matrix(Y0_g[match_c], nrow = nrow(Xc_s)), 1, stats::var)

  # Expectations by ingroup position
  expectations <- numeric(size)
  for (i in seq_len(size) - 1) {  # 0-indexed
    mask <- ingroup_labels == i
    arr <- var_t[mask] / (pg[mask] * q[mask]) +
      var_c[mask] / (pg[mask] * (1 - q[mask]))
    expectations[i + 1] <- mean(arr)
  }

  # Covariance of conditional treatment effects
  beta <- rowMeans(matrix(Y1_g[match_t], nrow = nrow(Xc_s))) -
    rowMeans(matrix(Y0_g[match_c], nrow = nrow(Xc_s)))
  beta_j <- lapply(seq_len(size) - 1, function(i) beta[ingroup_labels == i])
  min_len <- min(vapply(beta_j, length, integer(1)))
  beta_mat <- do.call(rbind, lapply(beta_j, function(b) b[seq_len(min_len)]))
  cov_mat <- stats::cov(t(beta_mat))

  Vg <- sum(expectations) / size^2 + sum(cov_mat) / size^2
  Vg
}


# ---- G encoding / decoding ----

#' Encode neighborhood structure G to integer
#' @keywords internal
encode_G_mat <- function(G, group_struct) {
  G <- as.matrix(G)
  weight <- c(1, cumprod(as.integer(group_struct[-length(group_struct)]) + 1))
  as.integer(G %*% weight)
}

#' Encode a single G vector
#' @keywords internal
encode_G_single <- function(g, group_struct) {
  weight <- c(1, cumprod(as.integer(group_struct[-length(group_struct)]) + 1))
  as.integer(sum(g * weight))
}

#' Decode integer to G vector
#' @keywords internal
decode_G <- function(G_encoded, group_struct) {
  weight <- c(1, cumprod(as.integer(group_struct[-length(group_struct)]) + 1))
  g <- integer(length(group_struct))
  remainder <- G_encoded
  for (i in length(weight):1) {
    g[i] <- remainder %/% weight[i]
    remainder <- remainder %% weight[i]
  }
  g
}


# ---- ClusterData builder ----

#' @keywords internal
cluster_data <- function(Y, Z, X, cluster_labels, group_labels, ingroup_labels,
                         cluster_feature, n_moments, categorical_Z) {
  Y <- as.numeric(Y)
  Z <- as.numeric(Z)
  X <- as.matrix(X)

  M <- length(unique(cluster_labels))
  units <- length(Y)

  # Split by group structure
  n_groups <- max(group_labels) + 1L
  n_clusters <- max(cluster_labels) + 1L

  # Build group_struct per cluster
  group_structs_mat <- matrix(0L, nrow = n_clusters, ncol = n_groups)
  for (i in seq_along(cluster_labels)) {
    group_structs_mat[cluster_labels[i] + 1L, group_labels[i] + 1L] <-
      group_structs_mat[cluster_labels[i] + 1L, group_labels[i] + 1L] + 1L
  }

  cluster_gs <- group_structs_mat[cluster_labels + 1L, , drop = FALSE]

  # Sort by group structure
  sort_order <- do.call(order, as.data.frame(cluster_gs))
  cluster_gs_sorted <- cluster_gs[sort_order, , drop = FALSE]

  data_by_group_struct <- list()
  idx1 <- 1
  first_group_struct <- NULL

  for (idx2 in seq_len(units)) {
    end_of_run <- (idx2 == units) ||
      !all(cluster_gs_sorted[idx2 + 1, ] == cluster_gs_sorted[idx1, ])

    if (end_of_run) {
      if (idx2 == units) {
        idx_range <- sort_order[idx1:idx2]
      } else {
        idx_range <- sort_order[idx1:(idx2)]
        # Check: only advance if not at end
      }

      gs <- cluster_gs_sorted[idx1, ]
      gs_name <- paste(gs, collapse = ",")

      if (is.null(first_group_struct)) first_group_struct <- gs

      sub_data <- get_final_tuple(
        Y[idx_range], Z[idx_range], X[idx_range, , drop = FALSE],
        cluster_labels[idx_range], group_labels[idx_range],
        ingroup_labels[idx_range],
        if (!is.null(cluster_feature)) cluster_feature[idx_range, , drop = FALSE] else NULL,
        gs, n_moments, categorical_Z)

      data_by_group_struct[[gs_name]] <- sub_data

      idx1 <- idx2 + 1
    }
  }

  list(
    data_by_group_struct = data_by_group_struct,
    M = M,
    first_group_struct = first_group_struct
  )
}


#' @keywords internal
get_final_tuple <- function(Y, Z, X, cluster_labels, group_labels,
                            ingroup_labels, cluster_feature, group_struct,
                            n_moments, categorical_Z) {
  # Sort by cluster label
  ord <- order(cluster_labels)
  Y <- Y[ord]
  Z <- Z[ord]
  X <- X[ord, , drop = FALSE]
  cluster_labels <- cluster_labels[ord]
  group_labels <- group_labels[ord]
  ingroup_labels <- ingroup_labels[ord]

  nunit_per_cluster <- sum(group_struct)
  k <- ncol(X)
  units <- length(Y)
  clusters <- units %/% nunit_per_cluster

  # Reshape X to 3D: clusters x units_per_cluster x k
  X_3d <- array(X, dim = c(nunit_per_cluster, clusters, k))
  X_3d <- aperm(X_3d, c(2, 1, 3))  # clusters x nunit x k

  X_aug <- X_3d

  if (n_moments > 0) {
    for (j in seq_along(group_struct)) {
      # Build mask for group j (0-indexed)
      group_j_mask <- matrix(group_labels, nrow = clusters, ncol = nunit_per_cluster, byrow = TRUE)
      group_j_mask <- (group_j_mask != (j - 1))  # TRUE = exclude

      # Mean of group j covariates per cluster
      Xj <- X_3d
      for (kk in seq_len(k)) {
        vals <- X_3d[, , kk]
        vals[group_j_mask] <- NA
        m1_k <- rowMeans(vals, na.rm = TRUE)  # clusters-length
        f1_k <- matrix(m1_k, nrow = clusters, ncol = nunit_per_cluster)
        X_aug <- abind_3d(X_aug, array(f1_k, dim = c(clusters, nunit_per_cluster, 1)))
      }

      # Higher moments
      if (n_moments > 1) {
        for (p in 2:n_moments) {
          for (kk in seq_len(k)) {
            vals <- X_3d[, , kk]
            vals[group_j_mask] <- NA
            m1_k <- rowMeans(vals, na.rm = TRUE)
            centered <- X_3d[, , kk] - matrix(m1_k, nrow = clusters, ncol = nunit_per_cluster)
            centered[group_j_mask] <- NA
            mp_k <- rowMeans(centered^p, na.rm = TRUE)
            fp_k <- matrix(mp_k, nrow = clusters, ncol = nunit_per_cluster)
            X_aug <- abind_3d(X_aug, array(fp_k, dim = c(clusters, nunit_per_cluster, 1)))
          }
        }
      }
    }
  }

  k_aug <- dim(X_aug)[3]
  # aperm to match Python's row-major reshape: units within each cluster must
  # be contiguous so that X_final rows align with Y, Z, G, and labels.
  X_final <- matrix(aperm(X_aug, c(2, 1, 3)), nrow = units, ncol = k_aug)

  # Compute G (number of treated neighbors per group, excluding self)
  ngroup <- length(group_struct)
  Z_onehot <- array(0, dim = c(clusters, nunit_per_cluster, ngroup))
  cl_idx <- cumsum(c(0, diff(cluster_labels)) > 0) + 1  # 1-indexed cluster index
  unit_in_cluster <- rep(seq_len(nunit_per_cluster), clusters)

  for (i in seq_along(Y)) {
    Z_onehot[cl_idx[i], unit_in_cluster[i], group_labels[i] + 1] <- Z[i]
  }

  G_plus <- array(NA, dim = c(clusters, nunit_per_cluster, ngroup))
  for (g in seq_len(ngroup)) {
    total_g <- rowSums(matrix(Z_onehot[, , g], nrow = clusters))
    G_plus[, , g] <- matrix(total_g, nrow = clusters, ncol = nunit_per_cluster)
  }
  G_stack <- G_plus - Z_onehot
  G <- matrix(aperm(G_stack, c(2, 1, 3)), nrow = units, ncol = ngroup)

  if (categorical_Z) {
    G <- matrix(as.integer(G), nrow = units, ncol = ngroup)
    Z <- as.integer(Z)
  }

  list(
    Y = Y, Z = Z, G = G, X = X_final,
    cluster_labels = cluster_labels,
    group_labels = group_labels,
    ingroup_labels = ingroup_labels,
    group_struct = group_struct
  )
}


#' Bind 3D arrays along third dimension
#' @keywords internal
abind_3d <- function(a, b) {
  d1 <- dim(a)
  d2 <- dim(b)
  result <- array(0, dim = c(d1[1], d1[2], d1[3] + d2[3]))
  result[, , seq_len(d1[3])] <- a
  result[, , d1[3] + seq_len(d2[3])] <- b
  result
}
