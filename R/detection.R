#' Calculate Feature Detection Rates
#'
#' Calculates feature-level detection rates using user-defined non-detect
#' values, feature-specific limits of detection (LODs), or both. The function
#' can optionally remove features whose detection rates are below a specified
#' threshold.
#'
#' @param data A data frame or matrix containing samples in rows and variables
#'   in columns. Matrix input is converted to a data frame.
#'
#' @param ID_col A single non-empty character string identifying the sample ID
#'   column. This argument must be provided explicitly.
#'
#' @param feature_cols An optional character vector containing the names of
#'   feature columns. When `NULL`, all columns except `ID_col` are treated as
#'   feature columns. If `ID_col` is included, it is removed with a warning.
#'
#' @param non_detect_values Optional values that should be counted as
#'   non-detects, such as `0`, `"0.00"`, `"ND"`, or `"<LOD"`.
#'   Numeric values are matched by value, so `0`, `"0"`, and `"0.00"` are
#'   treated the same. To count missing R values as non-detects, include
#'   an actual `NA`. If `NA` is not included, missing R values are skipped
#'   when calculating the detection rate.
#'
#' @param lod_table An optional data frame containing columns named `feature`
#'   and `lod`. The `feature` column identifies feature names and the `lod`
#'   column provides finite feature-specific LOD values. Negative LOD values
#'   are converted to zero with a warning. Additional columns are permitted
#'   and ignored.
#'
#' @param drop A single logical value indicating whether features should be
#'   removed according to their detection rates. When `FALSE`, detection rates
#'   are calculated without filtering. The default is `FALSE`.
#'
#' @param threshold A single numeric value between 0 and 1 used when
#'   `drop = TRUE`. Features with detection rates below this value, or with
#'   detection rates that cannot be calculated, are removed. Features with
#'   detection rates equal to the threshold are retained. The default is
#'   `0.5`.
#'
#' @details
#' At least one of `non_detect_values` or `lod_table` must be provided.
#'
#' When both are provided, a value is classified as detected only when it is
#' not a user-defined non-detect and is greater than or equal to its
#' feature-specific LOD.
#'
#' Actual `NA` values not included in `non_detect_values` are treated as
#' ordinary missing observations and excluded from the detection-rate
#' denominator.
#'
#' Negative feature measurements are converted to zero with a warning.
#' Negative LOD values are converted to zero with a warning.
#'
#' Duplicated input column names are reduced to their first occurrence.
#' Duplicated sample IDs are reduced to their first row. Both operations
#' produce warnings.
#'
#' @return A named list containing:
#'
#' \itemize{
#'   \item `data`: The cleaned data, optionally excluding dropped features.
#'   \item `detection_rates`: Feature-level LODs, evaluable counts, detected
#'     counts, detection rates, and filtering actions.
#'   \item `summary`: Numbers and proportion of retained and removed features.
#'   \item `data_issue`: Recorded data-quality messages.
#'   \item `kept_features`: Names of retained features.
#'   \item `dropped_features`: Names of removed features.
#' }
#'
#' @examples
#' example_data <- data.frame(
#'   sample_id = paste0("S", 1:4),
#'   feature_1 = c(0, 1.2, 2.1, NA),
#'   feature_2 = c(0, 0, 3.4, 4.2)
#' )
#'
#' detection(
#'   data = example_data,
#'   ID_col = "sample_id",
#'   non_detect_values = 0
#' )
#'
#' detection(
#'   data = example_data,
#'   ID_col = "sample_id",
#'   non_detect_values = 0,
#'   drop = TRUE,
#'   threshold = 0.5
#' )
#'
#' @export
#'
detection <- function(data,
                      ID_col,
                      feature_cols = NULL,
                      non_detect_values = NULL,
                      lod_table = NULL,
                      drop = FALSE,
                      # default: only detect
                      threshold = 0.5 # only works on drop = TRUE
){
  # Step 1. Check the main feature input data===================
  ## 1.1 Check the input data type==============================
  if (is.matrix(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or matrix.")
  }

  ### 1.2. Check duplicated column==========================
  duplicated_data_names <- unique(names(data)[duplicated(names(data))])
  #### 1.2.2.1 multiple id columns, retain only the first one.============================================================

  if (length(duplicated_data_names) > 0) {
    warning(
      "Duplicated column names were found in 'data': ",
      paste(duplicated_data_names, collapse = ", "),
      ". Only the first occurrence of each duplicated column will be retained."
    )

    data <- data[, !duplicated(names(data)), drop = FALSE]
  }

  ## 1.3 Check ID column========================================
  ### 1.3.1 Check ID column name type===========================
  if (missing(ID_col) ||
      !is.character(ID_col) ||
      length(ID_col) != 1 || is.na(ID_col) || trimws(ID_col) == "") {
    stop("'ID_col' must be provided as one non-empty column name.")
  }

  ### 1.3.2 Check ID column exists==============================
  if (!ID_col %in% names(data)) {
    stop(
      "The ID column '",
      ID_col,
      "' was not found in 'data'. Please check the name of the ID column."
    )
  }
  ### 1.3.3 Check missing ID values=============================
  if (anyNA(data[[ID_col]])) {
    stop("The ID column '", ID_col, "' contains missing values. ")
  }

  #### 1.3.3.2 Check duplicated IDs ==============================
  if (anyDuplicated(data[[ID_col]]) > 0) {
    duplicated_ids <- unique(data[[ID_col]][duplicated(data[[ID_col]])])
    warning(
      "Duplicated IDs were found in the ID column '",
      ID_col,
      "'. Only the first row for each duplicated ID will be retained. Duplicated IDs: ",
      paste(duplicated_ids, collapse = ", "),
      "."
    )

    # 'duplicated()' marks the second and later occurrence of each ID.
    data <- data[!duplicated(data[[ID_col]]), , drop = FALSE]
  }

  # Step 2. Drop options========================================
  ## 2.1 drop  = TRUE/FALSE; ===================================
  if (!is.logical(drop) || length(drop) != 1 || is.na(drop)) {
    stop(
      "'drop' must be either TRUE or FALSE.",
      " Note: TRUE means removing features below the detection-rate threshold. FALSE means calculating detection rates without removing any features. The default is FALSE."
    )
  }

  ## 2.2 check detection-rate threshold=========================
  if (drop) {
    # only when drop = TRUE
    if (!is.numeric(threshold) || length(threshold) != 1 ||
        is.na(threshold) || !is.finite(threshold) ||
        threshold < 0 || threshold > 1) {
      stop("'threshold' must be one real number between 0 and 1 when 'drop = TRUE'.")
    }
  }

  # Step 3. Identify feature columns============================
  if (is.null(feature_cols)) {
    ## 3.1 identify features using default======================
    # Default: all columns except the ID column as features
    feature_cols <- names(data)[names(data) != ID_col]

    if (length(feature_cols) == 0) {
      stop("No feature columns could be identified.")
    }
  }

  ## 3.2 User-provided features names===========================
  else {
    ### 3.2.1 Basic check=======================================
    if (!is.character(feature_cols) ||
        length(feature_cols) == 0 ||
        anyNA(feature_cols) || any(trimws(feature_cols) == "")) {
      stop(
        "Note: 'feature_cols' must be a non-empty character vector of column names without missing or blank values."
      )
    }

    # Remove the ID column if it is accidentally included in 'feature_cols'.
    if (ID_col %in% feature_cols) {
      warning("The ID column '",
              ID_col,
              "' was removed from 'feature_cols'.")

      feature_cols <- setdiff(feature_cols, ID_col)
    }

    ### 3.2.2 Remove duplicated requested feature names=========
    if (anyDuplicated(feature_cols) > 0) {
      duplicated_cols <- unique(feature_cols[duplicated(feature_cols)])
      warning(
        "Duplicated names in 'feature_cols' were removed: ",
        paste(duplicated_cols, collapse = ", "),
        "."
      )

      feature_cols <- unique(feature_cols)
    }

    # Removing the ID may leave no feature columns.
    if (length(feature_cols) == 0) {
      stop("No feature columns remain in 'feature_cols'.")
    }
  }

  ## 3.4 Every selected feature must exist in the input data.======
  missing_features <- setdiff(feature_cols, names(data))
  if (length(missing_features) > 0) {
    stop(
      "The following 'feature_cols' were not found in `data`: ",
      paste(missing_features, collapse = ", "),
      ". Please check the feature names"
    )
  }

  # Step 4. Prepare non-detect definitions======================
  ## 4.1 Require at least one: non_detect_values; lod_table.
  if ((is.null(non_detect_values) ||
       length(non_detect_values) == 0) && is.null(lod_table)) {
    stop("At least one of 'non_detect_values' or 'lod_table' must be provided.")
  }

  ## 4.2 User provided non_detect_values ========
  if (is.null(non_detect_values) ||
      length(non_detect_values) == 0) {
    ### 4.2.1 User only provide non_detect values.
    non_detect_includes_na <- FALSE
    non_detect_strings <- character(0)
    non_detect_numbers <- numeric(0)
  } else {
    ### 4.2.2 check 'non_detect_values' input type=======
    # make sure not a list, data frame, matrix, or array.
    if (!is.atomic(non_detect_values) ||
        !is.null(dim(non_detect_values))) {
      stop("'non_detect_values' must be an atomic vector.")
    }

    ### 4.2.3 handle R NA sepparately=============================
    # non_detect_values = c("NA")?
    non_detect_includes_na <- anyNA(non_detect_values)

    ### 4.2.4 store non-NA codes as unique character values.======
    # e.g.: 0; ND; <LOD
    non_detect_raw <- non_detect_values[!is.na(non_detect_values)]
    non_detect_strings <- unique(as.character(non_detect_raw))

    suppressWarnings(non_detect_numbers <- as.numeric(non_detect_strings))
    non_detect_numbers <- unique(non_detect_numbers[!is.na(non_detect_numbers)])
  }

  has_user_non_detect <- length(non_detect_strings) > 0 ||
    length(non_detect_numbers) > 0 ||
    non_detect_includes_na

  ## 4.3 Inspect actual missing ================================
  # create a empty character vector to store data-quality messages for inclusion in the returned object.
  data_issue_messages <- character(0)

  ### 4.3.1 Count actual R NA values in each selected feature=====
  actual_na_counts <- vapply(feature_cols, function(feature) {
    sum(is.na(data[[feature]]))
  }, numeric(1))

  n_actual_na_total <- sum(actual_na_counts)

  ### 4.3.2 Message if have actural NA==========================
  # to record how undeclared actual NA values will be handled
  if (n_actual_na_total > 0 && !non_detect_includes_na) {
    data_issue_messages <- c(
      data_issue_messages,
      paste0(
        "Actual NA values were detected in the feature data, but NA was not included in 'non_detect_values'. These values will remain missing and will be excluded from the detection-rate denominator."
      )
    )
  }

  ## 4.4 Inspect undeclared missing strings=====================
  ### 4.4.1 Define common missing strings=======================
  possible_missing_strings <- c("",
                                "NA",
                                "na",
                                "N/A",
                                "n/a",
                                "NaN",
                                "nan",
                                "NULL",
                                "null",
                                ".",
                                "missing",
                                "Missing")

  ### 4.4.2 Exclude if it is already defined by users===========
  possible_missing_strings_to_check <- setdiff(possible_missing_strings, non_detect_strings)
  # only remind possible_missing_strings - non_detect_strings
  #
  ### 4.4.3 Search each selected feature for undeclared missing-like strings.===========================================
  potential_missing_string_list <- lapply(feature_cols, function(feature) {
    raw_x <- data[[feature]]
    # Convert factors through character to avoid using internal factor codes.
    if (is.factor(raw_x)) {
      raw_x <- as.character(raw_x)
    }
    raw_string <- as.character(raw_x)

    potential_index <- !is.na(raw_string) &
      raw_string %in% possible_missing_strings_to_check
    # Return NULL when this feature has no potential missing strings.
    if (!any(potential_index)) {
      return(NULL)
    }

    # Count each potential missing string in the current feature.
    counts <- as.data.frame(table(raw_string[potential_index]), stringsAsFactors = FALSE)
    names(counts) <- c("potential_missing_string", "n")

    data.frame(
      feature = feature,
      potential_missing_string = as.character(counts$potential_missing_string),
      n = as.integer(counts$n),
      stringsAsFactors = FALSE
    )
  })

  # Remove NULL results from features without potential missing strings.
  potential_missing_string_list <- Filter(Negate(is.null), potential_missing_string_list)

  if (length(potential_missing_string_list) > 0) {
    # Combine feature-specific results into one table.
    potential_missing_strings <- do.call(rbind, potential_missing_string_list)

    row.names(potential_missing_strings) <- NULL

    data_issue_messages <- c(
      data_issue_messages,
      paste0(
        "Potential string-coded missing values were found but were not listed in 'non_detect_values'. Their feature names and values are recorded in 'potential_missing_strings'. Unrecognized non-numeric values will cause an error during numeric conversion."
      )
    )
  }

  else {
    # Return an empty table with a stable structure when no strings are found.
    potential_missing_strings <- data.frame(
      feature = character(0),
      potential_missing_string = character(0),
      n = integer(0),
      stringsAsFactors = FALSE
    )
  }


  # Step 5. Identify non-detects and covert to numeric. ========
  cleaned_data <- data
  user_non_detect_matrix <- vector("list", length(feature_cols))

  names(user_non_detect_matrix) <- feature_cols

  ## 5.1 Stop when Step 4 found undeclared missing-like strings=========================================================
  # Stop before numeric conversion if Step 4 found undeclared
  # missing-like strings.
  if (nrow(potential_missing_strings) > 0) {
    problem_details <- paste0(
      potential_missing_strings$feature,
      " = '",
      potential_missing_strings$potential_missing_string,
      "' (n = ",
      potential_missing_strings$n,
      ")"
    )
    stop(
      "Potential string-coded missing values were found but were not listed in 'non_detect_values': ",
      paste(problem_details, collapse = "; "),
      ". Add them to 'non_detect_values' if they represent non-detects, or convert them to actual NA if they represent missing observations."
    )
  }
  # Record the number of negative values found in each feature.
  negative_value_counts <- stats::setNames(integer(length(feature_cols)), feature_cols)
  ## 5.2 Identify user-defined non-detect observations==========
  for (feature in feature_cols) {
    # the orginial feature values
    raw_x <- cleaned_data[[feature]]
    # Convert factors through character first.
    if (is.factor(raw_x)) {
      raw_x <- as.character(raw_x)
    }
    # stop if unsupported onjects
    if (!(is.numeric(raw_x) ||
          is.integer(raw_x) || is.character(raw_x))) {
      stop(
        "Feature column '",
        feature,
        "' must contain numeric, integer, factor, or character values that can be converted to numeric."
      )
    }

    raw_string <- as.character(raw_x)

    # Identify actual R NA values
    actual_na_index <- is.na(raw_x)

    ## 5.3 to numeric and detect unknown==========================
    # Convert the feature values to numeric: declared string codes such as "ND" will become NA, but their positions remain recorded in `user_non_detect`.
    suppressWarnings(numeric_x <- as.numeric(raw_string))

    # Match only user-defined string codes
    string_match <- rep(FALSE, length(raw_x))

    if (length(non_detect_strings) > 0) {
      string_match <- !is.na(raw_string) &
        raw_string %in% non_detect_strings
    }

    # Match numeric non-detect values by numeric equivalence.
    numeric_match <- rep(FALSE, length(raw_x))

    if (length(non_detect_numbers) > 0) {
      numeric_match <- !is.na(numeric_x) &
        numeric_x %in% non_detect_numbers
    }

    # NA = actual NA and NA provided by user
    na_match <- actual_na_index & non_detect_includes_na

    # combine all user-defined non-detect
    # user_non_detect = match user defined non detect or (actual NA and NA provided by user)
    user_non_detect <- string_match | numeric_match | na_match



    # Find values that: 1. are not actual R NA; 2. were not declared as non-detects; and 3. cannot be converted to numeric. which cannot be interpreted
    bad_values <- unique(raw_string[!actual_na_index &
                                      !user_non_detect & is.na(numeric_x)])

    # Stop because the function cannot determine whether these values represent missing observations, non-detects, or data-entry errors.
    if (length(bad_values) > 0) {
      stop(
        "Feature column '",
        feature,
        "' contains unrecognized non-numeric value(s): ",
        paste(bad_values, collapse = ", "),
        ". Add them to 'non_detect_values' if they represent non-detects, or convert them to actual NA if they represent missing observations."
      )
    }

    ## Negative measured values are set to zero with an explicit warning.
    negative_index <- !is.na(numeric_x) & numeric_x < 0
    negative_value_counts[[feature]] <- sum(negative_index)
    if (any(negative_index)) {
      numeric_x[negative_index] <- 0
      # A converted zero must follow the user's zero non-detect rule.
      if (0 %in% non_detect_numbers) {
        user_non_detect[negative_index] <- TRUE
      }
    }

    ## 5.4 Store the numeric feature values.
    cleaned_data[[feature]] <- numeric_x

    # Store the corresponding non-detect value
    user_non_detect_matrix[[feature]] <- user_non_detect
  }

  features_with_negatives <- names(negative_value_counts[negative_value_counts > 0L])
  if (length(features_with_negatives) > 0L) {
    negative_details <- paste0(features_with_negatives,
                               " (n=",
                               negative_value_counts[features_with_negatives],
                               ")")
    negative_message <- paste0(
      "Negative feature values were converted to zero, following the meeting's provisional rule: ",
      paste(negative_details, collapse = ", "),
      ". Review the input if negative instrument values are meaningful."
    )
    warning(negative_message)
    data_issue_messages <- c(data_issue_messages, negative_message)
  }


  # Step 6. Check LOD table.====================================
  ## 6.0 Record whether the user supplied an LOD table.=========
  lod_provided <- !is.null(lod_table)

  ## 6.1 when no LOD table is provided.=========================
  lod_lookup <- numeric(0)
  ## 6.2 when LOD table is provided.============================
  if (lod_provided) {
    ### 6.2.1 Check the LOD table data type=====================
    if (!is.data.frame(lod_table)) {
      stop("'lod_table' must be a data frame.")
    }

    ### 6.2.2 Check duplicated column names ======================
    if (anyDuplicated(names(lod_table))) {
      duplicated_lod_columns <- unique(names(lod_table)[duplicated(names(lod_table))])

      stop(
        "'lod_table' contains duplicated column name(s): ",
        paste(duplicated_lod_columns, collapse = ", "),
        ". Column names in 'lod_table' must be unique."
      )
    }


    ### 6.2.3 Check the names of the LOD table: must use the fixed column names (feature, lod)===============================
    required_lod_cols <- c("feature", "lod")
    missing_lod_cols <- setdiff(required_lod_cols, names(lod_table))

    if (length(missing_lod_cols) > 0) {
      stop(
        "'lod_table' must contain columns named 'feature' and 'lod'. Missing column(s): ",
        paste(missing_lod_cols, collapse = ", "),
        "."
      )
    }

    ## 6.3 Retain only the columns needed for LOD matching.
    lod_clean <- lod_table[, required_lod_cols, drop = FALSE]

    ## 6.4 Check lod values=====================================
    ### 6.4.1 Convert factors to character.=====================
    if (is.factor(lod_clean$feature)) {
      lod_clean$feature <- as.character(lod_clean$feature)
    }

    # Feature names must be character values.
    lod_clean$feature <- as.character(lod_clean$feature)

    # Reject missing or blank feature names.
    if (any(is.na(lod_clean$feature) |
            trimws(lod_clean$feature) == "")) {
      stop("'lod_table$feature' contains missing or blank feature names.")
    }
    # Convert factor to character.
    if (is.factor(lod_clean$lod)) {
      lod_clean$lod <- as.character(lod_clean$lod)
    }
    original_lod <- lod_clean$lod
    ### 6.4.2 Convert factor to numeric.======================
    suppressWarnings(numeric_lod <- as.numeric(original_lod))

    ### 6.4.3 Identify values that could not be converted to numeric.========================================================
    invalid_lod <- is.na(numeric_lod) & !is.na(original_lod)
    if (any(invalid_lod)) {
      invalid_lod_values <- unique(original_lod[invalid_lod])
      stop(
        "'lod_table$lod' contains non-numeric value(s): ",
        paste(invalid_lod_values, collapse = ", "),
        "."
      )
    }

    ### 6.4.4 Check lod value data type:========================
    # # LOD values must be finite, non-missing.. # 0611 meeting

    # LOD values must be finite.
    if (any(!is.finite(numeric_lod))) {
      stop("'lod_table$lod' must contain non-missing, finite numeric values.")
    }

    negative_lod_index <- numeric_lod < 0

    if (any(negative_lod_index)) {
      negative_lod_message <- paste0(
        "Negative LOD values were converted to zero for feature(s): ",
        paste(unique(lod_clean$feature[negative_lod_index]), collapse = ", "),
        "."
      )

      warning(negative_lod_message)
      data_issue_messages <- c(data_issue_messages, negative_lod_message)

      numeric_lod[negative_lod_index] <- 0
    }
    lod_clean$lod <- numeric_lod
    ## 6.5. Check duplicated LOD feature names==================
    duplicated_lod_features <- unique(lod_clean$feature[duplicated(lod_clean$feature)])

    if (length(duplicated_lod_features) > 0) {
      # Identify duplicated features that have conflicting LOD values.
      conflicting_lod_features <- duplicated_lod_features[vapply(duplicated_lod_features, function(feature) {
        lod_values <- unique(lod_clean$lod[lod_clean$feature == feature])
        length(lod_values) > 1L
      }, logical(1))]

      # Conflicting LOD values are ambiguous and cannot be resolved safely.
      if (length(conflicting_lod_features) > 0L) {
        stop(
          "Different LOD values were provided for the same feature(s): ",
          paste(conflicting_lod_features, collapse = ", "),
          ". Each feature must have only one unique LOD value."
        )
      }

      # Reaching this point means duplicated rows have identical LOD values.
      warning(
        "Identical duplicated LOD rows were removed for: ",
        paste(duplicated_lod_features, collapse = ", "),
        "."
      )

      # Retain the first occurrence of each identical duplicated LOD row.
      lod_clean <- lod_clean[!duplicated(lod_clean$feature), , drop = FALSE]
    }

    # Create a named lookup vector: feature name -> LOD value.
    lod_lookup <- stats::setNames(lod_clean$lod, lod_clean$feature)

    ## 6.6 Check which selected features do not have an LOD value.=========================================================
    missing_lod_features <- setdiff(feature_cols, names(lod_lookup))

    if (length(missing_lod_features) > 0) {
      if (has_user_non_detect) {
        # Features without LOD values can still be evaluated using the user-defined non-detect codes.
        warning(
          "The following features do not have LOD values: ",
          paste(missing_lod_features, collapse = ", "),
          ". They will be evaluated using only 'non_detect_values'."
        )
      }

      else {
        # Without either an LOD or user-defined non-detect codes, the function has no rule for calculating detection rates for these features.
        warning(
          "The following features do not have LOD values, and no ",
          "'non_detect_values' were supplied: ",
          paste(missing_lod_features, collapse = ", "),
          ". Their detection rates will be returned as NA."
        )
      }
    }

    ## 6.7 check LOD table coverage ============================
    extra_lod_features <- setdiff(names(lod_lookup), feature_cols)
    if (length(extra_lod_features) > 0) {
      message(
        "'lod_table' contains LOD values for unselected feature(s): ",
        paste(extra_lod_features, collapse = ", "),
        ". These rows will be ignored."
      )
    }
  }

  # Step 7.Calculate feature-level detection rates==============
  ## 7.1 initialize result =====================================
  detection_results <- vector("list", length(feature_cols))
  names(detection_results) <- feature_cols

  ## 7.2 Process each selected feature==========================
  for (feature in feature_cols) {
    x <- cleaned_data[[feature]]
    user_non_detect <- user_non_detect_matrix[[feature]]

    # identify evaluable observations
    # Evaluable observations include: 1. numeric observations; and 2. values explicitly declared as non-detects.
    # Actual NA values not declared as non-detects are excluded.
    evaluable_index <- !is.na(x) | user_non_detect

    ### 7.2.1 LOD table provided===============
    ### With a feature-specific LOD: detected = evaluable numeric value at or above the LOD and not declared as a user-defined non-detect.
    has_lod <- lod_provided &&
      feature %in% names(lod_lookup)
    if (has_lod) {
      lod_value <- unname(lod_lookup[[feature]])
      # Define detect
      detected_index <- evaluable_index &
        !user_non_detect &
        !is.na(x) &
        x >= lod_value
    }
    # Current rule: values equal to LOD are treated as detected.
    # This boundary rule still requires final confirmation. # 0511 meeting
    #
    ### 7.2.2 Without LOD table================
    ### Without a feature-specific LOD: detected = evaluable numeric value not declared in non_detect_values
    # only user-defined values are classified as non-detects.
    else if (has_user_non_detect) {
      lod_value <- NA_real_
      detected_index <- evaluable_index &
        !user_non_detect &
        !is.na(x)
    }
    else {
      lod_value <- NA_real_
      evaluable_index <- rep(FALSE, length(x))
      detected_index <- rep(FALSE, length(x))
    }

    ## 7.3 Calculate the detection rate for each features========================================================
    n_evaluable <- sum(evaluable_index)
    n_detected <- sum(detected_index)

    # Return NA when no observations can be evaluated.
    detection_rate <- if (n_evaluable == 0) {
      NA_real_
    } else {
      n_detected / n_evaluable
    }

    # Store the feature result
    detection_results[[feature]] <- data.frame(
      feature = feature,
      lod = lod_value,
      n_evaluable = n_evaluable,
      n_detected = n_detected,
      detection_rate = detection_rate,
      stringsAsFactors = FALSE
    )
  }
  # Combine the feature-level results
  detection_rates <- do.call(rbind, detection_results)
  row.names(detection_rates) <- NULL

  # Step 8. Optionally drop low-detection feature columns======
  ## 8.1 Apply filtering when we choose drop = TRUE============
  if (drop) {
    # Features with a detection rate below the threshold and without an evaluable detection rate are removed.
    remove_index <- is.na(detection_rates$detection_rate) |
      detection_rates$detection_rate < threshold

    dropped_features <- detection_rates$feature[remove_index]
    kept_features <- detection_rates$feature[!remove_index]

    # which detection rate could not be calculated.
    missing_rate_features <- detection_rates$feature[is.na(detection_rates$detection_rate)]
    if (length(missing_rate_features) > 0) {
      warning(
        "The following features had no evaluable detection rate and were removed: ",
        paste(missing_rate_features, collapse = ", "),
        "."
      )
    }

    # Record the action
    detection_rates$action <- ifelse(remove_index, "removed", "kept")

    # Remove only the selected low-detection features.
    # ID and other non-feature columns remain unchanged.
    output_data <- cleaned_data[, setdiff(names(cleaned_data), dropped_features), drop = FALSE]
  }

  else {
    ## 8.2 Default: only calculate detection rates ========
    dropped_features <- character(0)
    kept_features <- detection_rates$feature
    detection_rates$action <- "not_applied"
    output_data <- cleaned_data
  }

  # Summarize the filtering results.
  n_total_features <- length(feature_cols)
  n_kept_features <- length(kept_features)
  n_dropped_features <- length(dropped_features)

  summary <- data.frame(
    n_total_features = n_total_features,
    n_kept_features = n_kept_features,
    n_dropped_features = n_dropped_features,
    kept_rate = n_kept_features / n_total_features,
    # Only when dropping was requested.
    threshold = if (drop)
      threshold
    else
      NA_real_,
    stringsAsFactors = FALSE
  )

  # Step 9. Return results======================================
  # Display recorded data-quality messages.
  if (length(data_issue_messages) > 0) {
    message(paste(unique(data_issue_messages), collapse = "\n"))
  }

  list(
    data = output_data,
    detection_rates = detection_rates,
    summary = summary,
    data_issue = unique(data_issue_messages),
    kept_features = kept_features,
    dropped_features = dropped_features
  )
}
