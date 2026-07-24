#' Replace values below LOD and/or non-detects in exposure data
#'
#' This function replaces values below limits of detection (LOD),
#' non-detect values (e.g. 0, NA), and user-defined missing codes
#' using flexible replacement strategies.
#'
#' Negative values are not allowed in either the input data or LOD table.
#' This typically indicates the data were scaled or transformed before
#' LOD correction.
#'
#' @param df A data frame containing exposure/chemical columns.
#' @param lod_table Optional data frame with columns `feature` and `lod`.
#' @param replacement_method Replacement method:
#' \itemize{
#'   \item{"lod_sqrt2"}{LOD / sqrt(2) (Default)}
#'   \item{"lod"}{LOD}
#'   \item{"zero"}{0}
#'   \item{"user"}{user_value}
#' }
#' @param user_value Numeric replacement value if method = "user".
#' @param id_col Optional ID column to exclude.
#' @param nondetect_codes Values treated as non-detects (e.g. c(0, NA)).
#' @param keep_missing If TRUE, nondetect_codes are kept as NA.
#' @param replace_below_lod If TRUE, values below LOD are replaced.
#' @param missing_codes Optional extra missing indicators (different from nondetects).
#' @param missing_code_action "replace" (with replacement_method) or "leave_missing"
#'
#' @return A data frame with transformed values.
#'
#' @export
replace_below_lod <- function(
    df,
    lod_table = NULL,
    replacement_method = c("lod_sqrt2", "lod", "zero", "user"),
    user_value = NULL,
    id_col = NULL,
    nondetect_codes = NULL,
    keep_missing = FALSE,
    replace_below_lod = TRUE,
    missing_codes = NULL,
    missing_code_action = c("replace", "leave_missing")
) {

  message("Starting LOD replacement function")

  replacement_method <- match.arg(replacement_method)
  missing_code_action <- match.arg(missing_code_action)

  cols_to_process <- setdiff(names(df), id_col)

  # ------------------------------------------------------------
  # REQUIREMENTS
  # ------------------------------------------------------------

  if (is.null(lod_table) && is.null(nondetect_codes)) {
    stop("Must provide either lod_table or nondetect_codes.")
  }

  if (replacement_method == "user" && is.null(user_value)) {
    stop("user_value required when replacement_method = 'user'.")
  }

  # ------------------------------------------------------------
  # VALIDATE LOD TABLE
  # ------------------------------------------------------------

  if (!is.null(lod_table)) {

    required_cols <- c("feature", "lod")

    if (!all(required_cols %in% names(lod_table))) {
      stop("lod_table must contain columns: feature and lod")
    }

    if (anyDuplicated(lod_table$feature)) {

      dupes <- unique(
        lod_table$feature[duplicated(lod_table$feature)]
      )

      stop(
        paste0(
          "Duplicate feature names in lod_table: ",
          paste(dupes, collapse = ", ")
        )
      )
    }

    if (any(!is.numeric(lod_table$lod))) {
      stop("All LOD values must be numeric.")
    }

    if (any(is.na(lod_table$lod))) {
      stop("LOD table contains missing LOD values.")
    }
  }

  # ------------------------------------------------------------
  # CHECK OVERLAPPING CODES
  # ------------------------------------------------------------

  if (!is.null(nondetect_codes) &&
      !is.null(missing_codes)) {

    overlap <- intersect(
      as.character(nondetect_codes),
      as.character(missing_codes)
    )

    if (length(overlap) > 0) {

      stop(
        paste0(
          "The following values appear in both ",
          "nondetect_codes and missing_codes: ",
          paste(overlap, collapse = ", ")
        )
      )
    }
  }

  # ------------------------------------------------------------
  # NEGATIVE CHECKS
  # ------------------------------------------------------------

  if (!is.null(lod_table) &&
      any(lod_table$lod < 0, na.rm = TRUE)) {

    stop("Negative values detected in LOD table.")
  }

  for (col in cols_to_process) {

    tmp <- suppressWarnings(as.numeric(df[[col]]))

    if (any(tmp < 0, na.rm = TRUE)) {
      stop(paste("Negative values detected in column:", col))
    }
  }

  # ------------------------------------------------------------
  # MESSAGES
  # ------------------------------------------------------------

  if (!is.null(id_col)) {
    message("ID column detected: ", id_col, " (excluded)")
  }

  message("Processing ", length(cols_to_process), " columns...")
  message("Replacement method: ", replacement_method)

  if (!is.null(nondetect_codes)) {
    message(
      "Non-detect codes: ",
      paste(nondetect_codes, collapse = ", ")
    )
  }

  if (!is.null(missing_codes)) {
    message(
      "Missing codes: ",
      paste(missing_codes, collapse = ", ")
    )
  }

  message("Keep missing: ", keep_missing)
  message("Replace below LOD: ", replace_below_lod)

  # ------------------------------------------------------------
  # CORE lapply
  # ------------------------------------------------------------

  result <- lapply(cols_to_process, function(col) {

    x_raw <- df[[col]]
    x <- suppressWarnings(
      as.numeric(as.character(x_raw))
    )

    # --------------------------
    # nondetects
    # --------------------------

    nondetect_idx <- rep(FALSE, length(x))

    if (!is.null(nondetect_codes)) {

      nondetect_idx <- as.character(x_raw) %in%
        as.character(nondetect_codes)

      if (any(is.na(nondetect_codes))) {
        nondetect_idx <- nondetect_idx | is.na(x)
      }
    }

    # --------------------------
    # missing codes
    # --------------------------

    missing_idx <- rep(FALSE, length(x))

    if (!is.null(missing_codes)) {

      missing_idx <- as.character(x_raw) %in%
        as.character(missing_codes)
    }

    # --------------------------
    # --------------------------
    # valid observed values
    # --------------------------

    exclude_idx <- nondetect_idx

    # exclude missing codes regardless of action
    # because they are not true observed values

    if (!is.null(missing_codes)) {
      exclude_idx <- exclude_idx | missing_idx
    }

    valid_positive <- x[
      !is.na(x) &
        x > 0 &
        !exclude_idx
    ]

    if (length(valid_positive) == 0) {
      stop(
        paste(
          "No valid positive observed values in:",
          col
        )
      )
    }

    min_positive <- min(valid_positive)

    # --------------------------
    # LOD
    # --------------------------

    lod_value <- NULL

    if (!is.null(lod_table)) {

      lod_match <- lod_table$lod[
        lod_table$feature == col
      ]

      # duplicate entries are not allowed
      if (length(lod_match) > 1) {

        stop(
          paste0(
            "Multiple LOD values found for feature: ",
            col
          )
        )
      }

      # no match OR NA -> fallback
      if (length(lod_match) == 0 || is.na(lod_match)) {

        warning(
          paste0(
            "Missing LOD value for feature: ",
            col,
            ". Using minimum observed positive value instead."
          )
        )

        lod_value <- min_positive

      } else {

        lod_value <- lod_match
      }
    }

    # --------------------------
    # replacement value
    # --------------------------

    base <- if (!is.null(lod_value)) {
      lod_value
    } else {
      min_positive
    }

    replacement_value <- switch(
      replacement_method,
      lod_sqrt2 = base / sqrt(2),
      lod = base,
      zero = 0,
      user = user_value
    )

    # --------------------------
    # counts
    # --------------------------

    below_count <- 0
    nondetect_count <- 0
    missing_count <- 0

    # --------------------------
    # below LOD
    # --------------------------

    if (replace_below_lod &&
        !is.null(lod_value)) {

      below_idx <- !is.na(x) &
        x < lod_value &
        !nondetect_idx &
        !missing_idx

      below_count <- sum(below_idx)

      x[below_idx] <- replacement_value
    }

    # --------------------------
    # nondetects
    # --------------------------

    if (keep_missing) {

      x[nondetect_idx] <- NA

    } else {

      nondetect_count <- sum(nondetect_idx)

      x[nondetect_idx] <- replacement_value
    }

    # --------------------------
    # missing codes
    # --------------------------

    if (!is.null(missing_codes)) {

      if (missing_code_action == "replace") {

        missing_count <- sum(missing_idx)

        x[missing_idx] <- replacement_value

      } else {

        x[missing_idx] <- NA
      }
    }

    list(
      data = x,
      below = below_count,
      nondetect = nondetect_count,
      missing = missing_count
    )
  })

  # ------------------------------------------------------------
  # UNPACK RESULTS
  # ------------------------------------------------------------

  df[cols_to_process] <- lapply(result, `[[`, "data")

  below_total <- sum(
    sapply(result, `[[`, "below")
  )

  nondetect_total <- sum(
    sapply(result, `[[`, "nondetect")
  )

  missing_total <- sum(
    sapply(result, `[[`, "missing")
  )

  # ------------------------------------------------------------
  # FINAL MESSAGES
  # ------------------------------------------------------------

  message(
    "Replaced ",
    below_total,
    " values below LOD."
  )

  if (!keep_missing) {

    message(
      "Replaced ",
      nondetect_total,
      " non-detect values."
    )
  }

  if (!is.null(missing_codes) &&
      missing_code_action == "replace") {

    message(
      "Replaced ",
      missing_total,
      " missing-code values."
    )
  }

  message("Replacement complete ✔")

  return(df)
}
