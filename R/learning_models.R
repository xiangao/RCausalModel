#' OLS Model Wrapper
#'
#' Creates a model wrapper around \code{lm()} that conforms to the CausalModel
#' model protocol: a list with \code{fit(X, y)} and \code{predict(model, X)}
#' functions.
#'
#' @return A list with \code{fit} and \code{predict} functions.
#' @export
cm_ols <- function() {
  list(
    fit = function(X, y) {
      X <- as.matrix(X)
      cnames <- paste0("V", seq_len(ncol(X)))
      colnames(X) <- cnames
      df <- data.frame(y = y, X)
      model <- stats::lm(y ~ ., data = df)
      model$.cm_colnames <- cnames
      model
    },
    predict = function(model, X) {
      X <- as.matrix(X)
      df <- data.frame(X)
      colnames(df) <- model$.cm_colnames
      as.numeric(stats::predict(model, newdata = df))
    },
    insample_predict = function(model) {
      as.numeric(stats::fitted(model))
    }
  )
}

#' Logistic Regression Model Wrapper
#'
#' Creates a model wrapper around \code{glm(family = binomial)} for binary
#' classification / propensity score estimation.
#'
#' @return A list with \code{fit}, \code{predict_proba}, and
#'   \code{insample_proba} functions.
#' @export
cm_logistic <- function() {
  list(
    fit = function(X, y) {
      X <- as.matrix(X)
      cnames <- paste0("V", seq_len(ncol(X)))
      colnames(X) <- cnames
      df <- data.frame(y = y, X)
      model <- stats::glm(y ~ ., data = df, family = stats::binomial())
      model$.cm_colnames <- cnames
      model
    },
    predict_proba = function(model, X) {
      X <- as.matrix(X)
      df <- data.frame(X)
      colnames(df) <- model$.cm_colnames
      as.numeric(stats::predict(model, newdata = df, type = "response"))
    },
    insample_proba = function(model) {
      as.numeric(stats::fitted(model))
    }
  )
}

#' Multinomial Logistic Regression Model Wrapper
#'
#' Creates a model wrapper around \code{nnet::multinom()} for multiclass
#' classification (e.g., neighborhood treatment propensity).
#'
#' @return A list with \code{fit}, \code{predict_proba}, and
#'   \code{insample_proba} functions.
#' @export
cm_multi_logistic <- function() {
  list(
    fit = function(X, y) {
      X <- as.matrix(X)
      cnames <- paste0("V", seq_len(ncol(X)))
      colnames(X) <- cnames
      df <- data.frame(y = as.factor(y), X)
      model <- nnet::multinom(y ~ ., data = df, trace = FALSE)
      model$.cm_colnames <- cnames
      model$classes <- sort(unique(y))
      model
    },
    predict_proba = function(model, X) {
      X <- as.matrix(X)
      df <- data.frame(X)
      colnames(df) <- model$.cm_colnames
      probs <- stats::predict(model, newdata = df, type = "probs")
      # If only 2 classes, multinom returns a vector; make it a matrix
      if (is.null(dim(probs))) {
        probs <- cbind(1 - probs, probs)
      }
      probs
    },
    insample_proba = function(model) {
      probs <- stats::fitted(model)
      if (is.null(dim(probs))) {
        probs <- cbind(1 - probs, probs)
      }
      probs
    }
  )
}

#' Random Forest Regressor Wrapper
#'
#' Creates a model wrapper around \code{randomForest::randomForest()} for
#' regression. Requires the \pkg{randomForest} package.
#'
#' @param ... Additional arguments passed to \code{randomForest::randomForest()}.
#' @return A list with \code{fit} and \code{predict} functions.
#' @export
cm_random_forest_regressor <- function(...) {
  rf_args <- list(...)
  list(
    fit = function(X, y) {
      if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("Package 'randomForest' is required. Install it with install.packages('randomForest').")
      }
      X <- as.matrix(X)
      args <- c(list(x = X, y = y), rf_args)
      model <- do.call(randomForest::randomForest, args)
      model$X_train <- X
      model
    },
    predict = function(model, X) {
      X <- as.matrix(X)
      as.numeric(stats::predict(model, newdata = X))
    },
    insample_predict = function(model) {
      as.numeric(stats::predict(model, newdata = model$X_train))
    }
  )
}

#' Random Forest Classifier Wrapper
#'
#' Creates a model wrapper around \code{randomForest::randomForest()} for
#' binary classification / propensity score estimation. Requires the
#' \pkg{randomForest} package.
#'
#' @param ... Additional arguments passed to \code{randomForest::randomForest()}.
#' @return A list with \code{fit}, \code{predict_proba}, and
#'   \code{insample_proba} functions.
#' @export
cm_random_forest_classifier <- function(...) {
  rf_args <- list(...)
  list(
    fit = function(X, y) {
      if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("Package 'randomForest' is required. Install it with install.packages('randomForest').")
      }
      X <- as.matrix(X)
      y_fac <- as.factor(y)
      args <- c(list(x = X, y = y_fac), rf_args)
      model <- do.call(randomForest::randomForest, args)
      model$X_train <- X
      model
    },
    predict_proba = function(model, X) {
      X <- as.matrix(X)
      probs <- stats::predict(model, newdata = X, type = "prob")
      as.numeric(probs[, 2])
    },
    insample_proba = function(model) {
      probs <- stats::predict(model, newdata = model$X_train, type = "prob")
      as.numeric(probs[, 2])
    }
  )
}
