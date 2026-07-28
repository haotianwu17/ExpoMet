#' Reformat Mass Spectrometry Feature Data
#'
#' Reformats mass spectrometry feature data from the common laboratory format,
#' in which features are stored in rows and samples are stored in columns, into
#' a more compatible format with samples in rows and features in columns.
#'
#' Feature-level metadata such as mass-to-charge ratio (`m/z`) and retention
#' time are preserved in a separate feature metadata table and linked to the
#' reformatted data using standardized feature IDs.
#'
#' @param data A data frame containing mass spectrometry feature data. Features
#'   should be stored in rows and samples in columns. The data frame must also
#'   contain columns identifying feature-level mass-to-charge ratio (`m/z`) and
#'   retention time.
#'
#' @param mz_col An optional single character string identifying the column
#'   containing mass-to-charge ratio (`m/z`) values. When `NULL`, the function
#'   attempts to identify the column automatically using common names such as
#'   `"mz"`, `"m/z"`, `"m.z"`, `"m_z"`, `"mass_to_charge"`, and
#'   `"mass.to.charge"`. If a unique column cannot be identified, an error is
#'   returned. The default is `NULL`.
#'
#' @param time_col An optional single character string identifying the column
#'   containing retention-time values. When `NULL`, the function attempts to
#'   identify the column automatically using common names such as `"time"`,
#'   `"rt"`, `"retention_time"`, and `"retention.time"`. If a unique column
#'   cannot be identified, an error is returned. The default is `NULL`.
#'
#' @param metadata_cols An optional character vector identifying additional
#'   feature-level metadata columns to retain, such as feature annotations,
#'   adducts, molecular formulas, ionization mode, or annotation confidence.
#'   These columns are stored in the returned feature metadata table and are
#'   not treated as sample columns. The `mz_col` and `time_col` columns are
#'   retained automatically and do not need to be included. The default is
#'   `NULL`.
#'
#' @param sample_cols An optional character vector identifying columns
#'   containing sample measurements. When `NULL`, all columns other than
#'   `mz_col`, `time_col`, and `metadata_cols` are treated as sample columns.
#'   The default is `NULL`.
#'
#' @param feature_prefix A single non-empty character string used as the prefix
#'   for newly generated feature IDs. Feature IDs are generated sequentially,
#'   for example `"F000001"`, `"F000002"`, and `"F000003"`. The default is
#'   `"F"`.
#'
#' @param ID_col A single non-empty character string specifying the name of the
#'   sample ID column in the reformatted data. Sample IDs are obtained from the
#'   original sample column names. The default is `"sample_id"`.
#'
#' @details
#' Mass spectrometry data are commonly exported with one feature per row,
#' feature characteristics such as `m/z` and retention time in separate
#' columns, and individual samples represented as additional columns. This
#' orientation is convenient for laboratory data processing but differs from
#' the conventional analysis format used in epidemiologic and statistical
#' workflows, where samples are represented in rows and variables are
#' represented in columns.
#'
#' `format_ms()` transposes the sample measurement matrix so that each sample
#' occupies one row and each mass spectrometry feature occupies one column.
#'
#' Standardized feature IDs are generated to provide stable identifiers for
#' features independent of their `m/z`, retention time, or annotation. The same
#' feature IDs are used as column names in the reformatted analytical data and
#' as identifiers in the returned feature metadata table.
#'
#' Feature-level information is therefore retained without requiring `m/z`,
#' retention time, or other metadata to be encoded directly into feature
#' names.
#'
#' When `sample_cols = NULL`, all columns not identified as feature metadata
#' are assumed to contain sample measurements. Users should explicitly provide
#' `metadata_cols` or `sample_cols` when the input data contain additional
#' non-sample columns.
#'
#' The function does not perform quality control, filtering, imputation,
#' normalization, transformation, or other preprocessing of feature
#' measurements. It only restructures the input data and preserves associated
#' feature metadata.
#'
#' @return A named list containing:
#'
#' \itemize{
#'   \item `data`: A data frame containing samples in rows and features in
#'     columns. The first column contains sample IDs using the name specified
#'     by `ID_col`.
#'   \item `feature_key`: A data frame containing the generated feature
#'     IDs, `m/z`, retention time, and any additional columns specified in
#'     `metadata_cols`.
#'   \item `sample_ids`: A character vector containing the sample IDs.
#'   \item `feature_ids`: A character vector containing the generated feature
#'     IDs.
#'   \item `summary`: A summary containing the number of samples, number of
#'     features, and retained feature metadata fields.
#' }
#'
#' @examples
#' example_data <- data.frame(
#'   mz = c(101.23, 342.18, 455.21),
#'   time = c(2.31, 5.72, 8.44),
#'   Sample_1 = c(100, 50, 210),
#'   Sample_2 = c(120, 45, 195),
#'   Sample_3 = c(98, 61, 220)
#' )
#'
#' formatted <- format_ms(example_data)
#'
#' formatted$data
#' formatted$feature_key
#'
#' example_annotated <- data.frame(
#'   mz = c(101.23, 342.18),
#'   rt = c(2.31, 5.72),
#'   adduct = c("[M+H]+", "[M-H]-"),
#'   annotation = c("Unknown", "Metabolite A"),
#'   Sample_1 = c(100, 50),
#'   Sample_2 = c(120, 45)
#' )
#'
#' formatted <- format_ms(
#'   data = example_annotated,
#'   mz_col = "mz",
#'   time_col = "rt",
#'   metadata_cols = c("adduct", "annotation")
#' )
#'
#' @export



format_ms <- function(data,
                      mz_col = NULL,
                      time_col = NULL,
                      metadata_cols = NULL,
                      sample_cols = NULL,
                      feature_prefix = "F",
                      ID_col = "sample_id") {
  
  # 1. Check input
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  
  if (!is.null(mz_col)) {
    if (length(mz_col) != 1 || !is.character(mz_col) ||
        is.na(mz_col) || mz_col == "") {
      stop("'mz_col' must be a single non-empty character string.")
    }
    if (!mz_col %in% names(data)) {
      stop("The column specified in 'mz_col' was not found in 'data'.")
    }
  }
  
  if (!is.null(time_col)) {
    if (length(time_col) != 1 || !is.character(time_col) ||
        is.na(time_col) || time_col == "") {
      stop("'time_col' must be a single non-empty character string.")
    }
    if (!time_col %in% names(data)) {
      stop("The column specified in 'time_col' was not found in 'data'.")
    }
  }
  
  if (!is.null(metadata_cols)) {
    if (!is.character(metadata_cols)) {
      stop("'metadata_cols' must be a character vector of column names.")
    }
    missing_metadata <- setdiff(metadata_cols, names(data))
    if (length(missing_metadata) > 0) {
      stop(
        "The following 'metadata_cols' were not found in 'data': ",
        paste(missing_metadata, collapse = ", ")
      )
    }
    metadata_cols <- unique(metadata_cols)
  }
  
  if (!is.null(sample_cols)) {
    if (!is.character(sample_cols)) {
      stop("'sample_cols' must be a character vector of column names.")
    }
    missing_samples <- setdiff(sample_cols, names(data))
    if (length(missing_samples) > 0) {
      stop(
        "The following 'sample_cols' were not found in 'data': ",
        paste(missing_samples, collapse = ", ")
      )
    }
    sample_cols <- unique(sample_cols)
  }
  
  if (length(feature_prefix) != 1 || !is.character(feature_prefix) ||
      is.na(feature_prefix) || feature_prefix == "") {
    stop("'feature_prefix' must be a single non-empty character string.")
  }
  
  if (length(ID_col) != 1 || !is.character(ID_col) ||
      is.na(ID_col) || ID_col == "") {
    stop("'ID_col' must be a single non-empty character string.")
  }
  
  

  # 2. Identify m/z column
  if (is.null(mz_col)) {
    mz_names <- c(
      "mz", "m/z", "m.z", "m_z",
      "mass_to_charge", "mass.to.charge",
      "mass-to-charge"
    )
    
    mz_match <- which(tolower(names(data)) %in% tolower(mz_names))
    
    if (length(mz_match) == 1) {
      mz_col <- names(data)[mz_match]
    } else if (length(mz_match) == 0) {
      stop(
        "Could not identify the m/z column automatically. ",
        "Specify it using 'mz_col'."
      )
    } else {
      stop(
        "Multiple possible m/z columns were identified: ",
        paste(names(data)[mz_match], collapse = ", "),
        ". Specify the correct column using 'mz_col'."
      )
    }
  }
  
  

  # 3. Identify retention-time column
  if (is.null(time_col)) {
    time_names <- c(
      "time", "rt",
      "retention_time", "retention.time",
      "retention time", "retention-time"
    )
    
    time_match <- which(tolower(names(data)) %in% tolower(time_names))
    
    if (length(time_match) == 1) {
      time_col <- names(data)[time_match]
    } else if (length(time_match) == 0) {
      stop(
        "Could not identify the retention-time column automatically. ",
        "Specify it using 'time_col'."
      )
    } else {
      stop(
        "Multiple possible retention-time columns were identified: ",
        paste(names(data)[time_match], collapse = ", "),
        ". Specify the correct column using 'time_col'."
      )
    }
  }
  
  

  # 4. Identify feature metadata and sample columns
  metadata_cols <- setdiff(metadata_cols, c(mz_col, time_col))
  
  feature_key_cols <- c(
    mz_col,
    time_col,
    metadata_cols
  )
  
  if (is.null(sample_cols)) {
    sample_cols <- setdiff(
      names(data),
      feature_key_cols
    )
  } else {
    overlap <- intersect(
      sample_cols,
      feature_key_cols
    )
    
    if (length(overlap) > 0) {
      stop(
        "The following columns were specified as both sample and ",
        "feature metadata columns: ",
        paste(overlap, collapse = ", ")
      )
    }
  }
  
  if (length(sample_cols) == 0) {
    stop("No sample columns were identified.")
  }
  
  non_numeric_samples <- sample_cols[
    !vapply(data[sample_cols], is.numeric, logical(1))
  ]
  
  if (length(non_numeric_samples) > 0) {
    stop(
      "All sample measurement columns must be numeric. ",
      "The following columns are not numeric: ",
      paste(non_numeric_samples, collapse = ", ")
    )
  }
  
  if (anyDuplicated(sample_cols)) {
    stop("Sample column names must be unique.")
  }
  
  
  # 5. Create standardized feature IDs
  n_features <- nrow(data)
  
  if (n_features == 0) {
    stop("'data' contains no feature rows.")
  }
  
  feature_ids <- paste0(
    feature_prefix,
    sprintf("%06d", seq_len(n_features))
  )
  

  # 6. Create feature metadata
  feature_key <- data.frame(
    feature_id = feature_ids,
    data[feature_key_cols],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(feature_key)[
    names(feature_key) == mz_col
  ] <- "mz"
  
  names(feature_key)[
    names(feature_key) == time_col
  ] <- "time"
  

  # 7. Transpose feature matrix
  measurement_matrix <- as.matrix(
    data[sample_cols]
  )
  
  formatted_matrix <- t(
    measurement_matrix
  )
  
  colnames(formatted_matrix) <- feature_ids
  

  # 8. Create analysis-ready data
  formatted_data <- data.frame(
    sample_id = sample_cols,
    formatted_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(formatted_data)[1] <- ID_col
  rownames(formatted_data) <- NULL
  
  # 9. Create summary
  cat(
    "format_ms() completed\n",
    "  Samples: ", length(sample_cols), "\n",
    "  Features: ", n_features, "\n",
    "  Feature metadata columns: ", length(feature_key_cols), "\n",
    sep = ""
  )

  # 10. Return results
  output <- list(
    data = formatted_data,
    feature_key = feature_key,
    sample_ids = sample_cols,
    feature_ids = feature_ids
  )
  
  return(output)
}

