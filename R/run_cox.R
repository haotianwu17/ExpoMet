#' Cox Proportional Hazard Model Analysis
#'
#' `run_cox` runs Cox proportional hazard regression models across all
#' exposure variables in the \code{omic_features} data frame. Models are fit using \code{survival::coxph}.
#' This function expects two dataframes: one data frame containing the omic features
#' data, specifying an ID column and omic columns, and one phenotype data
#' frame, containing the same ID column, a time-to-event column, an event
#' indicator column, and any covariates. Multiple testing correction is applied
#' using the Benjamini-Hochberg False Discovery Rate (FDR) method to identify
#' significant exposure-outcome associations.
#'
#' @param formula A two-sided formula specifying the full model, including the
#' survival outcome on the left-hand side using \code{Surv()}, e.g.
#' \code{Surv(time_to_event, event) ~ omic_features + age}. Variables from \code{pheno}
#' are referenced directly by column name. The keyword \code{omic_features} is used as
#' a placeholder for the exposure variables in the \code{omic_features} data frame,
#' which are automatically looped over. Interactions between omic features and phenotype
#' variables can be specified using \code{*} e.g. \code{Surv(time, event) ~ omic_features * covariate}.
#' When specifying an interaction, the main exposure of interest should always
#' be listed first. Only the interaction term coefficient is returned in the
#' results, not the main effect.
#' @param pheno A data frame containing the time-to-event column, event
#' indicator column, covariates, and optionally a \code{batch} column for
#' clustering. Must contain the ID column specified in \code{id_col}.
#' @param omic_features A data frame containing the exposure variables. Must
#' contain the ID column specified in \code{id_col}. Each column (excluding
#' the ID column) is treated as a separate exposure variable.
#' @param id_col A character string specifying the name of the ID column
#' present in both \code{pheno} and \code{omic_features} used to match participants
#' across the two data frames.
#' @param ties An optional character string specifying how to handle tied event
#' times. One of \code{"efron"} (default and recommended), \code{"breslow"},
#' or \code{"exact"}. See \code{?survival::coxph} for details. The exact method
#' is computationally expensive and not recommended for large datasets.
#' @param fdr.thres A numeric value specifying the Benjamini-Hochberg False
#' Discovery Rate (FDR) significance threshold.
#' Defaults to 0.05.
#' @param weights An optional vector or character string specifying the
#' regression weights. Defaults to NULL (no weights).
#' @param output_file An optional character string specifying a file path to
#' save the results table as a CSV file. Defaults to NULL (no file saved).
#'
#' @return A named list containing:
#' \describe{
#'   \item{results}{A data frame with the following columns:
#'   \itemize{
#'     \item \code{variable} - Name of the exposure variable from the
#'     \code{omic_features} data frame used in this model
#'     \item \code{term} - Name of main coefficient of interest extracted from
#'     the model
#'     \item \code{estimate} - Log hazard ratio (coefficient) from the fitted model
#'     \item \code{hazard_ratio} - Exponentiated estimate (hazard ratio)
#'     \item \code{std.error} - Standard error of the estimate
#'     \item \code{lcl} - Lower bound of confidence interval (log scale)
#'     \item \code{ucl} - Upper bound of confidence interval (log scale)
#'     \item \code{p.value} - Nominal p-value from the model
#'     \item \code{fdr} - Benjamini-Hochberg FDR adjusted p-value
#'     \item \code{significant} - Logical indicating FDR significance
#'   }}
#'   \item{plot}{A volcano plot displaying -log10(nominal p-value) on the
#'   y-axis against the log hazard ratio on the x-axis, with a Bonferroni
#'   threshold line indicated by a dashed horizontal line}
#' }
#'
# @importFrom survival coxph Surv cluster  # commented out as it will give error in building the package
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example without batch
#'
#' output_cox1 <- run_cox(
#'   Surv(time_to_event, event) ~ omic_features * ndi + age,
#'   pheno_df,
#'   omic_features,
#'   id_col      = "exp_id",
#'   output_file = "example_output/exwas_cox_results.csv"
#' )
#' output_cox1$results
#'
#' # Example with batch
#'
#' output_cox2 <- run_cox(
#'   Surv(time_to_event, event) ~ omic_features * ndi + age + cluster(batch),
#'   pheno_df,
#'   omic_features,
#'   id_col      = "exp_id",
#'   output_file = "example_output/exwas_cox_results.csv"
#' )
#'
#' }
#'
#'
#'


run_cox <- function(formula,
                    pheno,
                    omic_features,
                    id_col,
                    ties        = "efron",
                    fdr.thres   = 0.05,
                    weights     = NULL,
                    output_file = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(
      "Package 'survival' is required to use run_cox(). ",
      "Install it with install.packages('survival').",
      call. = FALSE
    )
  }

  # Step 1: Validate the inputs and merge data frames
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
  }
  if (length(formula) != 3) {
    stop("'formula' must be two-sided, with the survival outcome specified ",
         "on the left-hand side, e.g. Surv(time_to_event, event) ~ omic_features + age")
  }
  lhs_str <- deparse(formula[[2]])
  if (!grepl("^Surv\\(", trimws(lhs_str))) {
    stop("The left-hand side of 'formula' must be a Surv() call, e.g. ",
         "Surv(time_to_event, event) ~ omic_features + age. Do not prefix Surv() ",
         "with survival:: inside the formula.")
  }
  if (!is.data.frame(pheno)) {
    stop("'pheno' must be a data frame")
  }
  if (!is.data.frame(omic_features)) {
    stop("'omic_features' must be a data frame")
  }
  if (!id_col %in% colnames(pheno)) {
    stop("ID column '", id_col, "' not found in 'pheno'")
  }
  if (!id_col %in% colnames(omic_features)) {
    stop("ID column '", id_col, "' not found in 'omic_features'")
  }
  if (!ties %in% c("efron", "breslow", "exact")) {
    stop("'ties' must be one of 'efron', 'breslow', or 'exact'")
  }
  if (!is.numeric(fdr.thres) || fdr.thres <= 0 || fdr.thres >= 1) {
    stop("'fdr.thres' must be a number between 0 and 1")
  }

  # merge pheno and omic data by shared ID column
  data <- merge(pheno, omic_features, by = id_col)

  # extract exposure variable names (not ID column)
  exposures <- setdiff(colnames(omic_features), id_col)

  # check at least one exposure exists
  if (length(exposures) == 0) {
    stop("'omic_features' must contain at least one exposure column besides the ID")
  }

  # weights vector
  # validate and resolve weights
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1) {
      if (!weights %in% colnames(data)) {
        stop("Weights column '", weights, "' not found in 'pheno'. ")
      }
      w <- data[[weights]]
    } else if (is.numeric(weights)) {
      if (length(weights) != nrow(data)) {
        stop("Weights vector length (", length(weights), ") must match ",
             "the number of matched participants (", nrow(data), ")")
      }
      w <- weights
    } else {
      stop("'weights' must be either a column name (character string) ",
           "present in 'pheno', or a numeric vector of the same length ",
           "as the number of matched participants")
    }
  } else {
    w <- NULL
  }

  # check for participants lost during merge
  if (nrow(data) == 0) {
    stop("no matching participants found between 'pheno' and 'omic_features'
         — check that the ID column values match")
  }
  message(nrow(data), " participants matched across pheno and omic_features data frames")
  message(length(exposures), " exposure variables identified in omic_features")

  # collapse output to one line in case the formula is long
  formula_check   <- paste(deparse(formula), collapse = " ")
  has_interaction <- grepl("\\*", formula_check)

  errors   <- c()
  results  <- vector("list", length(exposures))
  progress <- txtProgressBar(min = 0, max = length(exposures), style = 3)

  for (i in seq_along(exposures)) {
    exp <- exposures[i]

    results[[i]] <- withCallingHandlers(
      tryCatch({

        # substitute the omic placeholder directly into the user's
        formula_str <- gsub("omic_features", exp, formula_check)

        model_formula <- as.formula(formula_str)

        model <- survival::coxph(model_formula, data = data, weights = w, ties = ties)

        coef_summary <- summary(model)$coefficients

        if (has_interaction) {
          coef_name <- rownames(coef_summary)[
            grep(paste0(exp, ":|:", exp), rownames(coef_summary))
          ]
          if (length(coef_name) == 0) stop("no interaction term found for exposure: ", exp)
        } else {
          coef_name <- exp
        }

        # se(coef) with no clustering, robust se with clustering
        se_col <- if ("robust se" %in% colnames(coef_summary)) "robust se" else "se(coef)"
        p_col  <- colnames(coef_summary)[grepl("^Pr", colnames(coef_summary))]

        result_list <- lapply(coef_name, function(coef_name) {
          est <- coef_summary[coef_name, "coef"]
          se  <- coef_summary[coef_name, se_col]

          data.frame(
            variable     = exp,
            term         = coef_name,
            estimate     = est,
            hazard_ratio = exp(est),
            std.error    = se,
            lcl          = est - 1.96 * se,
            ucl          = est + 1.96 * se,
            p.value      = coef_summary[coef_name, p_col]
          )
        })
        do.call(rbind, result_list)

      }, error = function(e) {
        errors[[exp]] <<- conditionMessage(e)
        NULL
      }),
      warning = function(w) {
        errors[[exp]] <<- paste0("Warning: ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    setTxtProgressBar(progress, i)
  }

  close(progress)

  if (length(errors) > 0) {
    message("\nThe following exposures failed or raised warnings and were skipped:")
    for (exp_name in names(errors)) message("  - ", exp_name, ": ", errors[[exp_name]])
  }

  results <- do.call(rbind, Filter(Negate(is.null), results))

  if (is.null(results) || nrow(results) == 0) {
    stop("No exposures produced a valid model. Check the errors reported above.")
  }

  results$fdr <- p.adjust(results$p.value, method = "fdr")
  results$significant <- results$fdr < fdr.thres
  results <- results[order(results$p.value), ]
  rownames(results) <- NULL

  log10_p <- -log10(results$p.value)
  bonferroni_thres <- 0.05 / length(exposures)
  n_significant <- sum(results$p.value < bonferroni_thres, na.rm = TRUE)

  if (n_significant == 0) {
    message("No exposures survived Bonferroni correction (threshold = ",
            signif(bonferroni_thres, 3), "). No significant hits will be displayed.")
  } else {
    message(n_significant, " exposure(s) survived Bonferroni correction.")
  }

  plot(
    x    = results$estimate,
    y    = log10_p,
    pch  = 16,
    xlab = "Log Hazard Ratio",
    ylab = expression(-log[10]*group("(", nominal~p-value, ")")),
    main = "Volcano Plot (Cox Prop.)"
  )
  abline(h = -log10(bonferroni_thres), col = "black", lty = 2, xpd = FALSE)
  volcano_plot <- recordPlot()

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    write.csv(results, file = output_file, row.names = FALSE)
    message("Results saved to: ", output_file)
  }

  return(list(results = results, plot = volcano_plot))
}

