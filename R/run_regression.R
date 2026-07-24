#' Univariate Regression Model Analysis
#'
#' `run_regression` runs regression models across all exposure variables in the
#' \code{omic_features} data frame. This function expects two dataframes; one
#' data frame containing the omic features data, specifying an ID column and omic feature
#' columns, and one phenotype data frame, containing the same ID column and
#' any covariates. Multiple testing correction is applied using the
#' Benjamini-Hochberg False Discovery Rate (FDR) method to identify significant
#' exposure-outcome associations.
#'
#'
#' @param formula A formula specifying the model. Variables from \code{pheno}
#' are referenced directly by column name. The keyword \code{omic_features} is used
#' as a placeholder for the exposure variables in \code{omic_features} data frame,
#' which are automatically looped over. Only the first predictor on the right-hand
#' side is extracted as the coefficient of interest. Remaining predictors are
#' treated as covariates to adjust for. A typical formula takes one of two forms:
#' \itemize{
#'   \item Omic features on right: \code{outcome ~ omic_features + covariates}
#'   \item Omic features on left: \code{omic_features ~ predictor + covariates}
#' }
#' Interactions between omic_features and phenotype variables can be specified
#' using \code{*} e.g. \code{outcome ~ omic_features * covariates}. When specifying an
#' interaction, the main coefficient of interest should always be listed
#' first in the formula (e.g. \code{outcome ~ omic_features * covariate}). Only the
#' interaction term coefficient is returned in the results, not the main effect.
#' @param pheno A data frame containing the outcome/predictor variable and
#' covariates. Must contain the ID column specified in \code{id}.
#' @param omic_features A data frame containing the exposure variables. Must
#' contain the ID column specified in \code{id}. Each column (excluding
#' the ID column) is treated as a separate exposure variable.
#' @param id_col A character string specifying the name of the ID column present
#' in both \code{pheno} and \code{omic_features} used to match participants
#' across the two data frames.
#' @param family A character string specifying the error distribution and link
#' function to be used in the function.
#' Use \code{"gaussian"} for an identity link, \code{"binomial"} for
#' a logit link, and \code{"poisson"} for a log link. Defaults to
#' \code{"gaussian"}.
#' @param fdr.thres A numeric value specifying the Benjamini-Hochberg False
#' Discovery Rate (FDR) significance threshold.
#' Defaults to 0.05.
#' @param weights An optional character string specifying the column name to be
#' used as regression weights. Defaults to NULL (no weights).
#' @param output_file An optional character string specifying a file path to
#' save the results table as a CSV file. Defaults to NULL (no file saved).
#'
#' @return A named list containing:
#' \describe{
#'   \item{results}{A data frame with the following columns:
#'   \itemize{
#'     \item \code{variable} - Name of the exposure variable from the
#'     \code{omic_features} data frame used in this regression model
#'     \item \code{term} - Name of main coefficient of interest extracted from
#'     the model
#'     \item \code{estimate} - Regression coefficient from the fitted model
#'     \item \code{std.error} - Standard error of the estimate
#'     \item \code{lcl} - Lower bound of confidence interval
#'     \item \code{ucl} - Upper bound of confidence interval
#'     \item \code{p.value} - Nominal p-value from the regression model
#'     \item \code{fdr} - Benjamini-Hochberg FDR adjusted p-value
#'     \item \code{significant} - Logical indicating FDR significance
#'   }}
#'   \item{plot}{A volcano plot displaying -log10(nominal p-value) on the
#'   y-axis against the effect estimate on the x-axis, with a Bonferroni
#'   threshold line indicated by a dashed horizontal line}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Displayed is an example of the function defined with an interaction term
#'
#' output2 <- exwas(
#'   outcome_case ~ omic_features * ndi + age,
#'   pheno_df,
#'   omic_features,
#'   id_col = "exp_id",
#'   family = "binomial",
#'   output_file = "example_output/exwas_interaction_results.csv"
#' )
#' output2$results
#' }
#'
#'
#'
exwas <- function(formula,
                    pheno,
                    omic_features,
                    id_col,
                    family      = "gaussian",
                    fdr.thres   = 0.05,
                    weights     = NULL,
                    output_file = NULL) {

  # Step 1: Validate the inputs and merge data frames
  # - check that formula is a glm formula object
  # - check that pheno and omic_features are data frames
  # - check that an ID column is found in pheno and omic_features
  # - check family is one of "binomial", "gaussian", or "poisson"
  # - check fdr.thres is between 0 and 1
  # - check that at least one exposure column exists in omic_features
  # - check for participants lost during merge step
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object")
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
  if (!family %in% c("gaussian", "binomial", "poisson")) {
    stop("'family' must be either 'gaussian', 'binomial', or 'poisson'")
  }
  if (!is.numeric(fdr.thres) || fdr.thres <= 0 || fdr.thres >= 1) {
    stop("'fdr.thres' must be a number between 0 and 1")
  }

  # merge pheno and omic_features by shared ID column
  data <- merge(pheno, omic_features, by = id_col)

  # extract exposure variable names (not ID column)
  exposures <- setdiff(colnames(omic_features), id_col)

  # check at least one exposure exists
  if (length(exposures) == 0) {
    stop("'omic_features' must contain at least one exposure column besides the ID")
  }

  # weights vector
  w <- if (!is.null(weights)) data[[weights]] else NULL

  # check for participants lost during merge
  if (nrow(data) == 0) {
    stop("no matching participants found between 'pheno' and 'omic_features'
         — check that the ID column values match")
  }
  message(nrow(data), " participants matched across pheno and omic_features data frames")
  message(length(exposures), " exposure variables identified in omic_features")

  # Step 2: Run regression models
  # - initialize a progress bar with style = 3 (to display percentage)
  # - for each exposure:
  #     - run glm() against the outcome
  #     - adjust for covariates if provided
  #     - extract estimate, std.error, ucl, lcl, and p.value from each model
  #     - update progress bar after each model
  # - close progress bar when all models are complete


  # detect which side omic_features is on
  formula_check <- deparse(formula)
  omic_features_on_left <- startsWith(trimws(formula_check), "omic_features")
  # detect whether there is an interaction term
  has_interaction <- grepl("\\*", deparse(formula))

  errors   <- c()    # store errors
  results  <- vector("list", length(exposures))
  progress <- txtProgressBar(min = 0, max = length(exposures), style = 3)

  for (i in seq_along(exposures)) {
    exp <- exposures[i]

    results[[i]] <- tryCatch({

      # build formula
      formula_str   <- gsub("omic_features", exp, formula_check)
      model_formula <- as.formula(formula_str)

      # fit model
      model <- glm(model_formula,
                   data    = data,
                   family  = family,
                   weights = w)

      # extract results
      coef_summary <- summary(model)$coefficients
      p_col        <- ifelse(family == "gaussian", "Pr(>|t|)", "Pr(>|z|)")

      if (omic_features_on_left) {
        first_predictor <- trimws(strsplit(deparse(formula[[3]]), "[+*]")[[1]][1])
        coef_name       <- first_predictor

        result_list <- lapply(coef_name, function(coef_name) {
          data.frame(
            variable  = exp,
            term      = coef_name,
            estimate  = coef_summary[coef_name, "Estimate"],
            std.error = coef_summary[coef_name, "Std. Error"],
            lcl       = coef_summary[coef_name, "Estimate"] - 1.96 * coef_summary[coef_name, "Std. Error"],
            ucl       = coef_summary[coef_name, "Estimate"] + 1.96 * coef_summary[coef_name, "Std. Error"],
            p.value   = coef_summary[coef_name, p_col]
          )
        })
        do.call(rbind, result_list)

      } else {
        if (has_interaction) {
          coef_name <- rownames(coef_summary)[
            grep(paste0(exp, ":|:", exp), rownames(coef_summary))
          ]
          if (length(coef_name) == 0) {
            stop("no interaction term found for exposure: ", exp)
          }
        } else {
          coef_name <- exp
        }

        result_list <- lapply(coef_name, function(coef_name) {
          data.frame(
            variable  = exp,
            term      = coef_name,
            estimate  = coef_summary[coef_name, "Estimate"],
            std.error = coef_summary[coef_name, "Std. Error"],
            lcl       = coef_summary[coef_name, "Estimate"] - 1.96 * coef_summary[coef_name, "Std. Error"],
            ucl       = coef_summary[coef_name, "Estimate"] + 1.96 * coef_summary[coef_name, "Std. Error"],
            p.value   = coef_summary[coef_name, p_col]
          )
        })
        do.call(rbind, result_list)
      }

    }, error = function(e) {
      # store error message for reporting
      errors[[exp]] <<- conditionMessage(e)
      NULL    # return NULL for this exposure so loop continues
    })

    setTxtProgressBar(progress, i)
  }

  close(progress)

  # report all errors after the loop completes
  if (length(errors) > 0) {
    message("\nThe following exposures failed and were skipped:")
    for (exp_name in names(errors)) {
      message("  - ", exp_name, ": ", errors[[exp_name]])
    }
  }

  # combine results into a single data frame
  results <- do.call(rbind, Filter(Negate(is.null), results))

  # Step 3: Apply FDR correction
  # - apply p.adjust(method = "fdr") to raw p-values
  # - flag exposures that meet the fdr.thres as significant
  # - sort results by fdr-adjusted p-value

  results$fdr <- p.adjust(results$p.value, method = "fdr")
  results$significant <- results$fdr < fdr.thres
  results <- results[order(results$p.value), ]
  rownames(results) <- NULL

  # Step 4: Generate volcano plot
  # x-axis: effect estimate (regression coefficient)
  # y-axis: -log10(nominal p-value)
  # - add horizontal line at bonferroni threshold value

  log10_p          <- -log10(results$p.value)

  bonferroni_thres <- 0.05 / length(exposures)

  n_significant <- sum(results$p.value < bonferroni_thres, na.rm = TRUE)
  if (n_significant == 0) {
    message("No exposures survived Bonferroni correction (threshold = ",
            signif(bonferroni_thres, 3), "). No significant hits will be ",
            "displayed on the volcano plot.")
  } else {
    message(n_significant, " exposure(s) survived Bonferroni correction.")
  }

  plot(
    x    = results$estimate,
    y    = log10_p,
    pch  = 16,
    xlab = "Effect Estimate (Regression Coefficient)",
    ylab = expression(-log[10]*group("(", nominal~p-value, ")")),
    main = "ExWAS Volcano Plot"
  )

  abline(h = -log10(bonferroni_thres), col = "black", lty = 2)

  volcano_plot <- recordPlot()


  # Step 5: Return named list
  # optionally save results to csv file
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    write.csv(results, file = output_file, row.names = FALSE)
    message("Results saved to: ", output_file)
  }
  return(list(results = results, plot = volcano_plot))

}


