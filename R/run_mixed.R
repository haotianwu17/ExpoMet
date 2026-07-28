#' Linear Mixed-Effects Model Analysis
#'
#' `run_mixed` runs linear mixed effects models across all exposure variables in the
#' \code{omic_feature} data frame, accounting for random intercepts and/or random slopes
#' (e.g. \code{(1 | batch)}). See lme4 for the appropiate string structure and syntax.
#' Models are fit using \code{lme4::lmer}
#' for a Gaussian family or \code{lme4::glmer} for binomial/Poisson families.
#' This function expects two dataframes; one data frame containing the omic features
#' data, specifying an ID column and omics columns, and one phenotype data frame,
#' containing the same ID column and any covariates. Multiple testing correction
#' is applied using the Benjamini-Hochberg False Discovery Rate (FDR) method to
#' identify significant exposure-outcome associations.
#'
#' @param formula A formula specifying the model, including fixed effects and random effects
#' in \code{lme4} syntax. Variables from \code{pheno}
#' are referenced directly by column name. The keyword \code{omic_features} is used
#' as a placeholder for the exposure variables in \code{omic_features} data frame,
#' which are automatically looped over. Only the first predictor on the right-hand
#' side is extracted as the coefficient of interest.
#' The formula should include at least one random effect term specified
#' using \code{lme4} syntax — e.g. \code{(1 | batch)} for a random
#' intercept by batch, or \code{(1 + age | batch)} for a random slope.
#' Random effects account for clustering or grouping structure in the data
#' (e.g. samples processed in the same batch sharing systematic measurement
#' variation). A typical formula takes the form:
#' \code{outcome ~ omic_features + covariates + (1 | grouping_variable)}.
#' Interactions between omic_features and phenotype variables can be
#' specified using \code{*}, e.g.
#' \code{outcome ~ omic_features * covariate + (1 | batch)}.
#' When specifying an interaction, the main exposure of interest should
#' always be listed first. Only the interaction term coefficient is returned
#' in the results, not the main effect

#' @param pheno A data frame containing the outcome/predictor variable, covariates,
#' and a \code{batch} column identifying cluster membership.
#' Must contain the ID column specified in \code{id_col}.
#' @param omic_features A data frame containing the exposure variables. Must
#' contain the ID column specified in \code{id_col}. Each column (excluding
#' the ID column) is treated as a separate exposure variable.
#' @param id_col A character string specifying the name of the ID column present
#' in both \code{pheno} and \code{omic_features} used to match participants
#' across the two data frames.
#' @param family A character string specifying the error distribution and link
#' function to be used in the function.
#' Use \code{"gaussian"} for an identity link, \code{"binomial"} for
#' a logit link, and \code{"Poisson"} for a log link. Defaults to
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
#' @importFrom lme4 lmer glmer
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # random intercept by batch — accounts for batch-level variation
#' run_mixed(
#'   outcome_case ~ omic_features + ndi + (1 | batch),
#'   pheno_df,
#'   omic_df,
#'   id_col = "exp_id",
#'   family = "binomial"
#' )
#'
#' # random slope — effect of age varies by batch
#' run_mixed(
#'   outcome_case ~ omic_features + age + (age | batch),
#'   pheno_df,
#'   omic_df,
#'   id_col = "exp_id",
#'   family = "binomial"
#' )
#'
#' }
#'
#'
#'
run_mixed <- function(formula,
                      pheno,
                      omic_features,
                      id_col,
                      family      = "gaussian",
                      fdr.thres   = 0.05,
                      weights     = NULL,
                      output_file = NULL) {

  # Step 1: Validate inputs
  if (!inherits(formula, "formula")) 
    stop("'formula' must be a formula object")
  
  if (!is.data.frame(pheno))    
    stop("'pheno' must be a data frame")
  
  if (!is.data.frame(omic_features)) 
    stop("'omic_features' must be a data frame")
  
  if (!id_col %in% colnames(pheno))    
    stop("ID column '", id_col, "' not found in 'pheno'")
  
  if (!id_col %in% colnames(omic_features)) 
    stop("ID column '", id_col, "' not found in 'omic_features'")
  
  if (!family %in% c("gaussian", "binomial", "poisson")) {
    stop("'family' must be either 'gaussian', 'binomial', or 'poisson'")
  }
  
  if (!is.numeric(fdr.thres) || fdr.thres <= 0 || fdr.thres >= 1) {
    stop("'fdr.thres' must be a number between 0 and 1")
  }

  formula_check <- paste(deparse(formula), collapse = " ")

  # merge pheno and omic_features by shared ID column
  data <- merge(pheno, omic_features, by = id_col)
  exposures <- setdiff(colnames(omic_features), id_col)
  if (length(exposures) == 0) stop("'omic_features' must contain at least one exposure column besides the ID")

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
  
  if (nrow(data) == 0) stop("no matching participants found between 'pheno' and 'omic_features'")
  message(nrow(data), " participants matched across pheno and omic_features data frames")
  message(length(exposures), " exposure variables identified in omic_features")

  omic_on_left <- startsWith(trimws(formula_check), "omic_features")
  has_interaction  <- grepl("\\*", formula_check)

  errors   <- c()
  results  <- vector("list", length(exposures))
  progress <- txtProgressBar(min = 0, max = length(exposures), style = 3)

  for (i in seq_along(exposures)) {
    exp <- exposures[i]

    results[[i]] <- withCallingHandlers(
      tryCatch({

        # build formula
        formula_str <- gsub("omic_features", exp, formula_check)
        model_formula     <- as.formula(formula_str)

        if (family == "gaussian") {
          model <- lme4::lmer(model_formula, data = data, weights = w)
        } else {
          model <- lme4::glmer(model_formula, data = data, family = family, weights = w)
        }

        coef_summary <- summary(model)$coefficients

        get_pvalue <- function(coef_name) {
          if (!coef_name %in% rownames(coef_summary)) {
            warning("coefficient '", coef_name, "' not found in model summary. Returning NA")
            return(NA_real_)
          }
          if (family == "gaussian") {
            t_val <- coef_summary[coef_name, "t value"]
            2 * pnorm(-abs(t_val))
          } else {
            coef_summary[coef_name, "Pr(>|z|)"]
          }
        }

        if (omic_on_left) {
          first_predictor <- trimws(strsplit(deparse(formula[[3]]), "[+*]")[[1]][1])
            coef_name <- rownames(coef_summary)[
              grep(paste0("^", first_predictor), rownames(coef_summary))
            ]
            if (length(coef_name) == 0) {
              stop("no coefficient found for predictor: ", first_predictor)
            }

          result_list <- lapply(coef_name, function(coef_name) {
            data.frame(
              variable  = exp,
              term      = coef_name,
              estimate  = coef_summary[coef_name, "Estimate"],
              std.error = coef_summary[coef_name, "Std. Error"],
              lcl       = coef_summary[coef_name, "Estimate"] - 1.96 * coef_summary[coef_name, "Std. Error"],
              ucl       = coef_summary[coef_name, "Estimate"] + 1.96 * coef_summary[coef_name, "Std. Error"],
              p.value   = get_pvalue(coef_name)
            )
          })
          do.call(rbind, result_list)

        } else {
          if (has_interaction) {
            coef_name <- rownames(coef_summary)[
              grep(paste0(exp, ":|:", exp), rownames(coef_summary))
            ]
            if (length(coef_name) == 0) stop("no interaction term found for exposure: ", exp)
          } else {
            coef_name <- rownames(coef_summary)[
              grep(paste0("^", exp), rownames(coef_summary))
            ]
            if (length(coef_name) == 0) {
              stop("no coefficient found for exposure: ", exp)
            }
          }

          result_list <- lapply(coef_name, function(coef_name) {
            data.frame(
              variable  = exp,
              term      = coef_name,
              estimate  = coef_summary[coef_name, "Estimate"],
              std.error = coef_summary[coef_name, "Std. Error"],
              lcl       = coef_summary[coef_name, "Estimate"] - 1.96 * coef_summary[coef_name, "Std. Error"],
              ucl       = coef_summary[coef_name, "Estimate"] + 1.96 * coef_summary[coef_name, "Std. Error"],
              p.value   = get_pvalue(coef_name)
            )
          })
          do.call(rbind, result_list)
        }

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
    xlab = "Effect Estimate (Regression Coefficient)",
    ylab = expression(-log[10]*group("(", nominal~p-value, ")")),
    main = "Volcano Plot (Mixed-Effects)"
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


