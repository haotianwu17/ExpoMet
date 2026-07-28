#' Format MWAS results for MetaboAnalyst Mummichog input
#'
#' This function formats MWAS (or similar metabolomics association results)
#' into a tab-delimited file compatible with MetaboAnalyst Mummichog pathway analysis.
#'
#' It supports:
#' \itemize{
#'   \item Optional grouping of results
#'   \item Effect size input (t-score)
#'   \item Optional retention time and ion mode annotation (e.g. pos or neg, C18, HILIC)
#' }
#'
#' @param exposure Character. Name of exposure or analysis (used for output file naming).
#'
#' @param input_df1 Data.frame containing MWAS results (e.g. C18 platform).
#'
#' @param input_df2 Optional second data.frame (e.g. HILIC platform). Exclude this if you only have one results data.frame
#'
#' @param annotation_df_list Data.frame OR list of data.frames containing feature annotations
#'   (e.g. c18_neg, c18_pos, hilic_neg, hilic_pos).
#'
#' @param output_dir Character. Local directory where formatted .txt files will be saved.
#'
#' @param feature_col Character. Feature identifier column name in the annotation_df file. Default is "feature_id".
#'
#' @param mz_col Character. Column name for mass-to-charge ratio. Default is "mz".
#'
#' @param rt_col Character or NULL. Column name for retention time.
#'
#' @param p_col Character. Column name for p-values. Default is "p.value".
#'
#' @param score_col Character. Column name for t-statistic or other effect size.
#'
#' @param group_col Character or NULL. Optional column used for stratified outputs.
#'
#' @param group_values Character vector or NULL. Optional values of group_col to iterate over.
#'
#' @param pos_neg_col Character or NULL. Column used to infer ionization mode.
#'
#' @return Invisibly returns TRUE. Writes formatted Mummichog input files to disk.
#'
#' @export

format_mummichog_input <- function(
    exposure,
    input_df1,
    input_df2 = NULL,
    annotation_df_list,
    output_dir,
    feature_col = "feature_id",
    mz_col = "mz",
    rt_col = NULL,
    p_col = "pvalue",
    score_col = "estimate",
    group_col = NULL,
    group_values = NULL,
    pos_neg_col = NULL
) {

  # -----------------------------
  # Validate inputs
  # -----------------------------
  if (!is.data.frame(input_df1)) {
    stop("input_df1 must be a data.frame")
  }

  if (!is.null(input_df2) && !is.data.frame(input_df2)) {
    stop("input_df2 must be a data.frame")
  }

  if (is.data.frame(annotation_df_list)) {
    dat_annot <- annotation_df_list
  } else if (is.list(annotation_df_list)) {

    if (!all(sapply(annotation_df_list, is.data.frame))) {
      stop("All elements of annotation_df_list must be data.frames")
    }

    dat_annot <- do.call(rbind, annotation_df_list)

  } else {
    stop("annotation_df_list must be a data.frame or list of data.frames")
  }

  if (is.null(score_col)) {
    stop("Must provide score_col.")
  }


  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # -----------------------------
  # Combine MWAS datasets
  # -----------------------------
  if (is.null(input_df2)) {
    dat_all <- input_df1
  } else {
    dat_all <- rbind(input_df1, input_df2)
  }

  # enforce ID column
  if (!(feature_col %in% colnames(dat_all))) {
    colnames(dat_all)[1] <- feature_col
  }

  # -----------------------------
  # Merge annotations
  # -----------------------------
  dat_all <- merge(dat_annot, dat_all, by = feature_col)

  # -----------------------------
  # Required columns
  # -----------------------------
  required_cols <- c(feature_col, mz_col, p_col, score_col)


  missing_cols <- setdiff(required_cols, colnames(dat_all))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # -----------------------------
  # Optional columns warning
  # -----------------------------
  optional_cols <- c(rt_col, pos_neg_col, group_col)
  optional_cols <- optional_cols[!is.null(optional_cols)]

  missing_optional <- setdiff(optional_cols, colnames(dat_all))

  if (length(missing_optional) > 0) {
    message(
      "Optional columns not found (ignored): ",
      paste(missing_optional, collapse = ", ")
    )
  }

  # -----------------------------
  # Ion mode
  # -----------------------------
  if (!is.null(pos_neg_col) &&
      pos_neg_col %in% colnames(dat_all)) {

    dat_all$mode <- ifelse(
      grepl("pos", dat_all[[pos_neg_col]], ignore.case = TRUE),
      "positive",
      "negative"
    )

  } else {
    dat_all$mode <- NA
  }

  # -----------------------------
  # Formatter
  # -----------------------------
  format_df <- function(df) {

    df <- df[order(df[[p_col]]), , drop = FALSE]

    score_name <- score_col

    out <- data.frame(
      m.z = df[[mz_col]],
      p.value = df[[p_col]],
      t.score = df[[score_name]]
    )

    if (!is.null(rt_col) && rt_col %in% colnames(df)) {
      out$r.t <- df[[rt_col]]
    }

    if ("mode" %in% colnames(df)) {
      out$mode <- df$mode
    }

    out
  }

  # -----------------------------
  # Write output
  # -----------------------------
  if (is.null(group_col)) {

    out <- format_df(dat_all)

    write.table(
      out,
      file.path(output_dir, paste0(exposure, "_mummichog.txt")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )

  } else {

    stopifnot(!is.null(group_values))

    for (g in group_values) {

      tmp <- subset(dat_all, dat_all[[group_col]] == g)
      out <- format_df(tmp)

      write.table(
        out,
        file.path(output_dir, paste0(exposure, "_", g, "_mummichog.txt")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  }

  invisible(TRUE)
}
