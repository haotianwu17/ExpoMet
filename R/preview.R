#' Preview of Data
#'
#'@description
#' `preview` can be used as a quick overview of common data issues that may
#' affect model fitting and is designed to be used prior to any subsequent
#' function in the ExpoMet package. This function checks and identifies
#' duplicate IDs, non-numeric columns, data missingness, negative numbers, exposures with zeroes,
#' exposures with zero variance, exposures with significantly skewed data, using
#' the  D'Agostino-Pearson test (FDR < 0.05), and whether the minimum value in the dataset is
#' overrepresented by more than five times what would be expected by chance.
#' The first 5 features are returned for each check, with the option to view the
#' entire dataset.
#'
#'
#' @param omic_features A dataframe with an id column and exposure columns.
#' @param id_col A character string specifying the name of the id column
#' in \code{omic_features}.


#' @return A printed summary containing:
#'
#' \describe{
#' \item{Duplicate IDs}{identifies any duplicate IDs present in the ID column.}
#' \item{Non-numeric columns}{the exposure columns that contain non-numeric values
#' and the percentage of columns in the dataset with non-numeric values.}
#'  \item{Missing data}{the exposure columns and rows that contain missing values (\code{NA}
#'     or \code{NaN}) and the percentage of missing values per column or row.}
#'  \item{Zeroes Present}{the exposure columns containing zeroes and the percentage
#'  of zeroes per column.}
#'  \item{Zero variance}{the exposure columns with zero variance across all samples.}
#'  \item{Skewness}{the exposure columns with skewness that significantly deviates from
#'  a normal distribution, evaluated using the D'Agostino-Pearson test.}
#'  \item{Minimum Inflation}{the exposure columns that have an overrepresented minimum value
#'  of more than five times what is expected by chance.}
#' }
#'
#' @examples
#' \dontrun{
#' preview(omic_features, id_col = "exp_id")
#'}
#'
#' @export

preview <- function(omic_features, id_col) {

  # check before running any computation
  if (!is.data.frame(omic_features)) {
    stop("'omic_features' must be a data frame")
  }
  if (!id_col %in% names(omic_features)) {
    stop("ID column '", id_col, "' not found in 'omic_features'. ")
  }
  if (nrow(omic_features) == 0) {
    stop("'omic_features' contains no rows")
  }
  exposure_cols <- setdiff(names(omic_features), id_col)
  if (length(exposure_cols) == 0) {
    stop("No exposure columns found in 'omic_features' after removing the ID column '",
         id_col, "'")
  }

  # define skewness test
  d_agostino_skew_test <- function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 8) return(c(skewness = NA, p_value = NA))

    xbar <- mean(x)
    m2 <- mean((x - xbar)^2)
    m3 <- mean((x - xbar)^3)
    g1 <- m3 / m2^(3/2)

    # D'Agostino transformation
    Y <- g1 * sqrt((n + 1) * (n + 3) / (6 * (n - 2)))
    beta2 <- 3 * (n^2 + 27*n - 70) * (n + 1) * (n + 3) /
      ((n - 2) * (n + 5) * (n + 7) * (n + 9))
    W2 <- -1 + sqrt(2 * (beta2 - 1))
    delta <- 1 / sqrt(log(sqrt(W2)))
    alpha <- sqrt(2 / (W2 - 1))
    Z <- delta * log(Y / alpha + sqrt((Y / alpha)^2 + 1))

    p_value <- 2 * pnorm(-abs(Z))

    c(skewness = g1, p_value = p_value)
  }


# get the exposure columns only (no id column)

exposure_df <- omic_features[exposure_cols]
numeric_df  <- exposure_df[sapply(exposure_df, is.numeric)]

# Step 1: identify duplicate IDs
id_values <- omic_features[[id_col]]
id_table  <- table(id_values)

duplicate_id_summary <- data.frame(
  id    = names(id_table[id_table > 1]),
  count = as.integer(id_table[id_table > 1])
)
duplicate_id_summary <- duplicate_id_summary[order(-duplicate_id_summary$count), ]
rownames(duplicate_id_summary) <- NULL

# Step 2: Identify Non-numeric Columns
is_non_numeric  <- !sapply(exposure_df, is.numeric)
non_numeric_summary <- data.frame(exposure = names(exposure_df)[is_non_numeric])
non_numeric_pct <- round(nrow(non_numeric_summary) / length(exposure_cols) * 100, 2)

# Step 3: Identify Missing Data (columns and rows)
# cols
col_missing_pct <- sapply(exposure_df, function(x)
  suppressWarnings(mean(is.na(x) | is.nan(as.numeric(x))))
)
col_missing_df <- data.frame(
  column      = names(col_missing_pct),
  col_missing = round(col_missing_pct * 100, 2)
)
col_missing_df <- col_missing_df[col_missing_df$col_missing > 0, ]
col_missing_df <- col_missing_df[order(-col_missing_df$col_missing), ]

# rows
row_missing_pct <- apply(exposure_df, 1, function(x)
  suppressWarnings(mean(is.na(x) | is.nan(as.numeric(x))))
)
row_missing_df <- data.frame(
  row         = omic_features[[id_col]],
  row_missing = round(row_missing_pct * 100, 2)
)
row_missing_df <- row_missing_df[row_missing_df$row_missing > 0, ]
row_missing_df <- row_missing_df[order(-row_missing_df$row_missing), ]

# combined into one summary table
n_col <- nrow(col_missing_df)
n_row <- nrow(row_missing_df)
n_max <- max(n_col, n_row)

missing_summary <- data.frame(
  column      = c(col_missing_df$column,      rep(NA, n_max - n_col)),
  col_missing = c(col_missing_df$col_missing, rep(NA, n_max - n_col)),
  row         = c(row_missing_df$row,         rep(NA, n_max - n_row)),
  row_missing = c(row_missing_df$row_missing, rep(NA, n_max - n_row))
)

# Step 4: Identify Negative Numbers
negative_pct <- sapply(numeric_df, function(x) mean(x < 0, na.rm = TRUE))
negative_summary <- data.frame(
  exposure = names(negative_pct),
  negative_pct = round(negative_pct * 100, 2),
  row.names = NULL
)

negative_summary <- negative_summary[negative_summary$negative_pct > 0, ]
negative_summary <- negative_summary[order(-negative_summary$negative_pct), ]
rownames(negative_summary) <- NULL

# Step 5: Identify Zeroes
zero_pct <- sapply(numeric_df, function(x) mean(x == 0, na.rm = TRUE))
zero_summary <- data.frame(
  exposure = names(zero_pct),
  zero_pct = round(zero_pct * 100, 2),
  row.names = NULL
)

zero_summary <- zero_summary[zero_summary$zero_pct > 0, ]
zero_summary <- zero_summary[order(-zero_summary$zero_pct), ]
rownames(zero_summary) <- NULL

# Step 6: Identify No Variance
variances <- sapply(numeric_df, var, na.rm = TRUE)
variance_summary <- data.frame(
  exposure = names(variances)[!is.na(variances) & variances == 0],
    row.names = NULL)

# Step 7: Skewness Indicator (with FDR correction)
skew_results <- sapply(numeric_df, d_agostino_skew_test)

skew_summary <- data.frame(
  exposure  = colnames(skew_results),
  skewness  = round(skew_results["skewness", ], 2),
  p_value   = signif(skew_results["p_value", ], 3),
  row.names = NULL
)

skew_summary <- skew_summary[!is.na(skew_summary$p_value), ]
skew_summary$fdr <- p.adjust(skew_summary$p_value, method = "fdr")

# filter on FDR-adjusted p-value
skew_summary <- skew_summary[skew_summary$fdr < 0.05, ]
skew_summary <- skew_summary[order(skew_summary$fdr), ]
rownames(skew_summary) <- NULL

# Step 8: Minimum Inflation
min_inflation <- sapply(numeric_df, function(x){
  x <- x[!is.na(x)]
  if (length(x) == 0) return(c(min_value = NA, observed_freq = NA, fold_over_expected = NA))

  min_val <- min(x)

  if (min_val == 0) return(c(min_value = min_val, observed_freq = NA, fold_over_expected = NA))

  nonNA_n       <- length(x)
  expected_freq <- 1 / nonNA_n

  observed_freq <- mean(x == min_val)
  fold_over_expected <- observed_freq / expected_freq

  c(min_value = min_val, observed_freq = observed_freq, fold_over_expected = fold_over_expected)
})

min_inflation_summary <- data.frame(
  exposure            = colnames(min_inflation),
  min_value           = min_inflation["min_value", ],
  observed_pct        = round(min_inflation["observed_freq", ] * 100, 2),
  fold_over_expected  = round(min_inflation["fold_over_expected", ], 1),
  row.names = NULL
)

min_inflation_summary <- min_inflation_summary[!is.na(min_inflation_summary$fold_over_expected), ]
min_inflation_summary <- min_inflation_summary[min_inflation_summary$fold_over_expected >= 5, ]
min_inflation_summary <- min_inflation_summary[order(-min_inflation_summary$fold_over_expected), ]
rownames(min_inflation_summary) <- NULL

# Print Results for User
cat("Quality Control Data Summary \n\n")
cat("Duplicate IDs\n")
if (nrow(duplicate_id_summary) == 0) {
  cat("No duplicate IDs detected.\n\n")
} else {
  cat(nrow(duplicate_id_summary), "duplicate ID(s) identified.\n")
  print(head(duplicate_id_summary, 5))
  if (nrow(duplicate_id_summary) > 5) {
    cat("See `duplicate_id` in the returned object for all results.\n")
  }
  cat("\n")
}

cat("Non-Numeric Columns\n")
if (nrow(non_numeric_summary) == 0) {
  cat("No non-numeric columns detected.\n\n")
} else {
  cat(nrow(non_numeric_summary), "non-numeric column(s) detected",
      paste0("(", non_numeric_pct, "% of all exposure columns):\n"))
  print(head(non_numeric_summary, 5))
  if (nrow(non_numeric_summary) > 5) {
    cat("See 'non_numeric' in the returned object for all results.\n")
  }
  cat("\n")
}

cat("Missing Data\n")
if (n_col == 0 && n_row == 0) {
  cat("No missing data detected.\n\n")
} else {
  cat(n_col, "column(s) and", n_row, "row(s) with missing data:\n")
  print(head(missing_summary, 5))
  if (nrow(missing_summary) > 5) {
    cat("See 'missing' in the returned object for all results.\n")
  }
  cat("\n")
}

cat("Negative Numbers\n")
if (nrow(negative_summary) == 0) {
  cat("No negative numbers detected.\n\n")
} else {
  cat(nrow(negative_summary), "exposure(s) with negative numbers:\n")
  print(head(negative_summary, 5))
  if (nrow(negative_summary) > 5) {
    cat("See 'negative' in the returned object for all results.\n")
  }
  cat("\n")
}

cat("Zeroes Present\n")
if (nrow(zero_summary) == 0) {
  cat("No zeroes detected.\n\n")
} else {
  cat(nrow(zero_summary), "exposure(s) with zeroes:\n")
  print(head(zero_summary, 5))
  if (nrow(zero_summary) > 5) {
    cat("See 'zero' in the returned object for all results.\n")
  }
  cat("\n")
}

cat("Zero Variance\n")
if (nrow(variance_summary) == 0) {
  cat("No zero variance exposures detected.\n\n")
} else {
  cat(nrow(variance_summary), "exposure(s) with zero variance:\n")
  print(head(variance_summary, 5))
  if (nrow(variance_summary) > 5) {
    cat("See 'variance' in the returned object for all results.\n")
  }
  cat("\n")
}

cat("Skewness\n")
if (nrow(skew_summary) == 0) {
  cat("No significantly skewed features detected (FDR < 0.05).\n\n")
} else {
  cat(nrow(skew_summary), "exposure(s) with significant skewness (FDR < 0.05):\n")
  print(head(skew_summary, 5))
  if (nrow(skew_summary) > 5) {
    cat("See 'skewness' in the returned object for all skewed features.\n")
  }
  cat("\n")
}

cat("Minimum Value Inflation\n")
if (nrow(min_inflation_summary) == 0) {
  cat("No overrepresented minimum values detected.\n\n")
} else {
  cat(nrow(min_inflation_summary), "exposure(s) with overrepresented minimum values:\n")
  print(head(min_inflation_summary, 5))
  if (nrow(min_inflation_summary) > 5) {
    cat("See 'min_inflation' in the returned object for all results.\n")
  }
  cat("\n")
}

invisible(list(
  duplicate_id  = duplicate_id_summary,
  non_numeric   = non_numeric_summary,
  missing       = missing_summary,
  negative      = negative_summary,
  zero          = zero_summary,
  variance      = variance_summary,
  skewness      = skew_summary,
  min_inflation = min_inflation_summary
))

}

