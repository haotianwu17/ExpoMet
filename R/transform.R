#' Transform Feature Columns
#'
#' Applies a selected transformation method to user-specified feature columns
#' while preserving the sample ID column and other non-feature columns.
#'
#' This function is intended to be used after detection and replacement steps,
#' and before downstream statistical analysis. It does not perform non-detect
#' detection, non-detect replacement, imputation, or outlier correction.
#'
#' @param data A data frame or matrix containing samples in rows and variables
#'   in columns. Matrix input is converted to a data frame.
#'
#' @param id_col A single non-empty character string identifying the sample ID
#'   column. This column is never transformed.
#'
#' @param feature_cols An optional character vector containing the names of
#'   feature columns to transform. When \code{NULL}, all columns except
#'   \code{id_col} are treated as feature columns. If \code{id_col} is included,
#'   it is removed with a warning. If the input data contain covariate columns,
#'   users should explicitly provide \code{feature_cols} to avoid transforming
#'   non-feature variables.
#'
#' @param method A character string specifying the transformation method.
#'   Options are:
#'   \itemize{
#'     \item \code{"log2"}: log2 transformation.
#'     \item \code{"zscore"}: z-score scaling.
#'     \item \code{"pareto"}: Pareto scaling.
#'     \item \code{"boxcox"}: Box-Cox transformation.
#'   }
#'
#' @details
#' Let \eqn{x_{ij}} denote the value of feature \eqn{j} for sample \eqn{i},
#' where \eqn{i = 1, \ldots, n} and \eqn{j = 1, \ldots, p}.
#'
#' For log2 transformation, each selected feature value is transformed as:
#'
#' \deqn{
#' x^{*}_{ij} = \log_2(x_{ij})
#' }
#'
#' This transformation requires \eqn{x_{ij} > 0}. If any selected feature
#' column contains values less than or equal to zero, the function stops and
#' reports all problematic feature columns. Users should run the upstream
#' detection and replacement functions before applying log2 transformation.
#'
#' For z-score scaling, each selected feature is centered by its sample mean
#' and scaled by its sample standard deviation:
#'
#' \deqn{
#' x^{*}_{ij} = \frac{x_{ij} - \bar{x}_j}{s_j}
#' }
#'
#' where
#'
#' \deqn{
#' \bar{x}_j = \frac{1}{n_j}\sum_{i \in O_j} x_{ij}
#' }
#'
#' and
#'
#' \deqn{
#' s_j =
#' \sqrt{
#' \frac{1}{n_j - 1}
#' \sum_{i \in O_j}
#' (x_{ij} - \bar{x}_j)^2
#' }
#' }
#'
#' Here, \eqn{O_j} is the set of non-missing observations for feature
#' \eqn{j}, and \eqn{n_j = |O_j|}.
#'
#' For Pareto scaling, each selected feature is centered by its sample mean
#' and scaled by the square root of its sample standard deviation:
#'
#' \deqn{
#' x^{*}_{ij} = \frac{x_{ij} - \bar{x}_j}{\sqrt{s_j}}
#' }
#'
#' For Box-Cox transformation, each selected feature is transformed using a
#' feature-specific parameter \eqn{\lambda_j}. When \eqn{\lambda_j \ne 0}:
#'
#' \deqn{
#' x^{*}_{ij} =
#' \frac{x_{ij}^{\lambda_j} - 1}{\lambda_j}
#' }
#'
#' When \eqn{\lambda_j = 0}:
#'
#' \deqn{
#' x^{*}_{ij} = \log(x_{ij})
#' }
#'
#' Box-Cox transformation requires \eqn{x_{ij} > 0}. If any selected feature
#' column contains values less than or equal to zero, the function stops and
#' reports all problematic feature columns. The parameter \eqn{\lambda_j} is
#' estimated separately for each feature by maximizing the Box-Cox
#' log-likelihood over the interval \eqn{[-2, 2]}.
#'
#' Columns with zero or missing standard deviation are not transformed by
#' z-score or Pareto scaling.
#'
#' @return A data frame with the selected feature columns transformed using
#'   the specified method. The sample ID column and other non-feature columns
#'   are preserved unchanged.
#'
#' @examples
#' example_data <- data.frame(
#'   sample_id = paste0("S", 1:4),
#'   feature_1 = c(1.2, 2.4, 3.1, 4.8),
#'   feature_2 = c(5.5, 6.2, 7.1, 8.3),
#'   age = c(30, 35, 40, 45)
#' )
#'
#' transform_features(
#'   data = example_data,
#'   id_col = "sample_id",
#'   feature_cols = c("feature_1", "feature_2"),
#'   method = "log2"
#' )
#'
#' transform_features(
#'   data = example_data,
#'   id_col = "sample_id",
#'   feature_cols = c("feature_1", "feature_2"),
#'   method = "zscore"
#' )
#'
#' transform_features(
#'   data = example_data,
#'   id_col = "sample_id",
#'   feature_cols = c("feature_1", "feature_2"),
#'   method = "pareto"
#' )
#'
#' transform_features(
#'   data = example_data,
#'   id_col = "sample_id",
#'   feature_cols = c("feature_1", "feature_2"),
#'   method = "boxcox"
#' )
#'
#' @export

transform_features <- function(data,
                               id_col,
                               feature_cols = NULL,
                               method = c("log2",
                                          "zscore",
                                          "pareto",
                                          "boxcox")) {
  method <- match.arg(method)

  if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or matrix.")
  }

  duplicated_data_names <- unique(names(data)[duplicated(names(data))])

  if (length(duplicated_data_names) > 0) {
    warning(
      "Duplicated column names were found in 'data': ",
      paste(duplicated_data_names, collapse = ", "),
      ". Only the first occurrence of each duplicated column will be retained."
    )

    data <- data[, !duplicated(names(data)), drop = FALSE]
  }

  if (missing(id_col) ||
      !is.character(id_col) ||
      length(id_col) != 1 ||
      is.na(id_col) ||
      trimws(id_col) == "") {
    stop("'id_col' must be provided as one non-empty column name.")
  }

  if (!id_col %in% names(data)) {
    stop(
      "The ID column '",
      id_col,
      "' was not found in 'data'. Please check the name of the ID column."
    )
  }

  if (anyNA(data[[id_col]])) {
    stop("The ID column '", id_col, "' contains missing values.")
  }

  if (anyDuplicated(data[[id_col]]) > 0) {
    duplicated_ids <- unique(data[[id_col]][duplicated(data[[id_col]])])

    warning(
      "Duplicated IDs were found in the ID column '",
      id_col,
      "'. Only the first row for each duplicated ID will be retained. Duplicated IDs: ",
      paste(duplicated_ids, collapse = ", "),
      "."
    )

    data <- data[!duplicated(data[[id_col]]), , drop = FALSE]
  }

  if (is.null(feature_cols)) {
    feature_cols <- names(data)[names(data) != id_col]

    if (length(feature_cols) == 0) {
      stop("No feature columns could be identified.")
    }
  } else {
    if (!is.character(feature_cols) ||
        length(feature_cols) == 0 ||
        anyNA(feature_cols) ||
        any(trimws(feature_cols) == "")) {
      stop("'feature_cols' must be a non-empty character vector of column names without missing or blank values.")
    }

    if (id_col %in% feature_cols) {
      warning("The ID column '", id_col, "' was removed from 'feature_cols'.")
      feature_cols <- setdiff(feature_cols, id_col)
    }

    if (anyDuplicated(feature_cols) > 0) {
      duplicated_cols <- unique(feature_cols[duplicated(feature_cols)])

      warning(
        "Duplicated names in 'feature_cols' were removed: ",
        paste(duplicated_cols, collapse = ", "),
        "."
      )

      feature_cols <- unique(feature_cols)
    }

    if (length(feature_cols) == 0) {
      stop("No feature columns remain in 'feature_cols'.")
    }
  }

  missing_features <- setdiff(feature_cols, names(data))

  if (length(missing_features) > 0) {
    stop(
      "The following 'feature_cols' were not found in 'data': ",
      paste(missing_features, collapse = ", "),
      ". Please check the feature names."
    )
  }

  non_numeric_cols <- feature_cols[
    !vapply(data[feature_cols], is.numeric, logical(1))
  ]

  if (length(non_numeric_cols) > 0) {
    stop(
      "All `feature_cols` must be numeric. Non-numeric columns: ",
      paste(non_numeric_cols, collapse = ", ")
    )
  }

  if (method %in% c("log2", "boxcox")) {
    non_positive_cols <- feature_cols[
      vapply(
        data[feature_cols],
        function(x) any(x <= 0, na.rm = TRUE),
        logical(1)
      )
    ]

    if (length(non_positive_cols) > 0) {
      stop(
        "The following feature columns contain values <= 0 and cannot be transformed using '",
        method,
        "': ",
        paste(non_positive_cols, collapse = ", "),
        ". Please run the detection and replacement functions before transformation."
      )
    }
  }

  if (method == "boxcox") {
    zero_variance_cols <- feature_cols[
      vapply(
        data[feature_cols],
        function(x) {
          x_non_missing <- x[!is.na(x)]
          length(unique(x_non_missing)) < 2
        },
        logical(1)
      )
    ]

    if (length(zero_variance_cols) > 0) {
      stop(
        "The following feature columns have zero or insufficient variance and cannot be transformed using 'boxcox': ",
        paste(zero_variance_cols, collapse = ", "),
        ". Box-Cox transformation requires each selected feature column to contain at least two distinct non-missing values."
      )
    }
  }

  transformed_data <- data

  for (col in feature_cols) {
    x <- data[[col]]

    if (method == "log2") {
      # Values <= 0 have already been checked before the loop.
      transformed_data[[col]] <- log2(x)
    }

    if (method == "zscore") {
      col_mean <- mean(x, na.rm = TRUE)
      col_sd <- stats::sd(x, na.rm = TRUE)

      if (is.na(col_sd) || col_sd == 0) {
        warning(
          "Column `", col,
          "` has zero or missing standard deviation. ",
          "Z-score scaling was not applied to this column."
        )
      } else {
        transformed_data[[col]] <- (x - col_mean) / col_sd
      }
    }

    if (method == "pareto") {
      col_mean <- mean(x, na.rm = TRUE)
      col_sd <- stats::sd(x, na.rm = TRUE)

      if (is.na(col_sd) || col_sd == 0) {
        warning(
          "Column `", col,
          "` has zero or missing standard deviation. ",
          "Pareto scaling was not applied to this column."
        )
      } else {
        transformed_data[[col]] <- (x - col_mean) / sqrt(col_sd)
      }
    }

    if (method == "boxcox") {
      # Values <= 0 have already been checked before the loop.
      x_non_missing <- x[!is.na(x)]

      boxcox_transform <- function(value, lambda) {
        if (abs(lambda) < 1e-8) {
          log(value)
        } else {
          (value^lambda - 1) / lambda
        }
      }

      boxcox_loglik <- function(lambda) {
        y <- boxcox_transform(x_non_missing, lambda)

        n <- length(y)
        variance <- mean((y - mean(y))^2)

        if (is.na(variance) || variance <= 0) {
          return(-Inf)
        }

        -n / 2 * log(variance) +
          (lambda - 1) * sum(log(x_non_missing))
      }

      optimal_lambda <- stats::optimize(
        f = function(lambda) -boxcox_loglik(lambda),
        interval = c(-2, 2)
      )$minimum

      transformed_data[[col]] <- boxcox_transform(x, optimal_lambda)
    }
  }

  transformed_data
}
