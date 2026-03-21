#' Create a Potential Outcome Data Container
#'
#' Validates and stores (Y, Z, X) data, splitting into treated/control groups.
#'
#' @param Y Numeric vector. Outcome variable.
#' @param Z Numeric vector. Treatment indicator (0/1 for binary, continuous allowed).
#' @param X Numeric matrix or vector. Covariates. Vectors are converted to
#'   single-column matrices.
#' @return A \code{po_data} object (S3 list).
#' @export
po_data <- function(Y, Z, X = NULL) {
  Y <- as.numeric(Y)
  Z <- as.numeric(Z)
  n <- length(Y)

  if (length(Z) != n) {
    stop("Y and Z must have the same length")
  }

  if (!is.null(X)) {
    if (is.vector(X) || is.factor(X)) {
      X <- matrix(as.numeric(X), ncol = 1)
    } else {
      X <- as.matrix(X)
    }
    if (nrow(X) != n) {
      stop("X must have the same number of rows as length of Y")
    }
  }

  idx_t <- which(Z == 1)
  idx_c <- which(Z == 0)

  structure(
    list(
      Y = Y,
      Z = Z,
      X = X,
      n = n,
      idx_t = idx_t,
      idx_c = idx_c,
      nt = length(idx_t),
      nc = length(idx_c),
      Yt = Y[idx_t],
      Yc = Y[idx_c],
      Xt = if (!is.null(X)) X[idx_t, , drop = FALSE] else NULL,
      Xc = if (!is.null(X)) X[idx_c, , drop = FALSE] else NULL
    ),
    class = "po_data"
  )
}
