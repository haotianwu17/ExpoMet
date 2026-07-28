#' Exploratory Data Analysis Summary
#'
#'@description
#' `exp_eda` is designed to provide a summary of common exploratory data
#' analysis steps for exposome-wide or metabolomic-wide association studies. This function identifies
#' and visualizes the distribution of outliers in feature intensities across samples, performs
#' Principal Component analysis across the omic features, optionally displays visual
#' differences by batch and tests their differences using ANOVA, and identifies correlated
#' features.
#'
#' @param omic_features A dataframe with an id column and omic features columns.
#' @param id_col A character string specifying the name of the id column
#' in \code{omic_features}.
#' @param rho_thresh A numeric value specifying the minimum Spearman correlation
#' threshold for two features to be grouped in the same correlation block. The value
#' must be between 0 and 1. The default is \code{0.5}.
#' @param min_block_size An integer specifying the minimum number of features required
#' for a correlation block to be created and reported. Blocks with fewer than this
#' threshold are discarded. Default is \code{5}.
#' @param batch An optional vector of batch assignments, one per row in
#'   \code{omic_features}. If provided, stacked boxplots of PC1 and PC2 colored by batch
#'   will be displayed. The PCA scatterplot will be colored by batch, and ANOVA
#'   will be used to test for significant batch effects on PC1 and PC2. Default is
#'   \code{NULL} (no batch coloring or testing).
#'
#' @return A combination of printed summaries and plots, including:
#'
#' \describe{
#'   \item{Forest Plot Relative Intensity Summary}{Individual boxplots showing the
#'     distribution of relative intensities across samples for those with outliers, displayed from minimum
#'     intensity to maximum intensity. Outliers (values beyond 3 standard deviations
#'     from the mean) are shown in unfilled black circles. All other values are
#'     reflected only in the box and whiskers, not as individual points.}
#'   \item{PCA Variance Table}{A printed table showing the proportion of
#'     variance and cumulative variance explained by the top 5 principal
#'     components, computed on complete cases across non-zero-variance
#'     exposure columns.}
#'     \item{PCA Top 10 Absolute Loadings Table}{A table showing the top 10
#'     absolute loadings for the top 5 principal components, identifying which
#'     individual features have the strongest overall influence on the component.
#'     The first five rows are shown in the console and the full table is available
#'     in the returned object.}
#'     \item{Batch Effect Assessment}{An assessment of significant batch effects using ANOVA
#'     (significant defined as p < 0.05).}
#'     \item{PC1 and PC2 Boxplots}{If a batch vector is provided, stacked boxplots
#'     showing the distribution of relative intensity by batch for PC1 and PC2 will
#'     be displayed.}
#'   \item{Scree Plot}{A bar chart showing the percentage of variance
#'     explained by each of the first 15 principal components, if applicable.}
#'   \item{PCA Scatterplot}{A scatterplot of the first two principal
#'     components, with axis labels showing the percent variance each
#'     explains. If a batch vector is provided, the scatterplot is colored by batch.}
#'   \item{Correlation Blocks Table}{A printed table of the top correlated features, calculated using Spearman rank correlation,
#'      identified through hierarchical clustering. The first five blocks
#'      are shown in the console and the full table is available in the returned object.}
#'
#' }
#'
#' @examples
#' set.seed(1)
#'
#' omic_features <- data.frame(
#'   exp_id = 1:20,
#'   A = rnorm(20),
#'   B = rnorm(20),
#'   C = rnorm(20)
#' )
#'
#' batch <- rep(c("Batch1", "Batch2"), each = 10)
#' 
#' exp_eda(omic_features,
#'         id_col         = "exp_id",
#'         rho_thresh      = 0.5,
#'         min_block_size = 6,
#'         batch          = batch
#'         )
#'
#' @export
#'

exp_eda <- function(omic_features, id_col, rho_thresh = 0.5,
                    min_block_size = 5, batch = NULL) {

  exposure_cols <- setdiff(names(omic_features), id_col)
  numeric_df <- omic_features[exposure_cols]
  numeric_df <- numeric_df[sapply(numeric_df, is.numeric)]

  if (!is.null(batch)) {
    if (!is.vector(batch) && !is.factor(batch)) {
      stop("'batch' must be a vector of batch assignments, ",
           "one entry per row in 'omic_features'")
    }
    if (length(batch) != nrow(omic_features)) {
      stop("'batch' vector length (", length(batch), ") must match ",
           "the number of rows in 'omic_features' (", nrow(omic_features), ")")
    }
    batch <- as.factor(batch)

    if (nlevels(batch) < 2) {
      warning("'batch' contains only one unique value. ANOVA and coloring ",
              "by batch require at least 2 batches. Batch will be ignored.")
      batch <- NULL
    }
  }

  if (!is.numeric(rho_thresh) || rho_thresh <= 0 || rho_thresh >= 1) {
    stop("'rho_thresh' must be a number between 0 and 1")
  }
  if (!is.numeric(min_block_size) || min_block_size < 2) {
    stop("'min_block_size' must be an integer of at least 2")
  }

# Step 1: Forest Plots - Relative Intensity Summary Per Experiment

  # filtered to only samples with at least one outlier
  sample_ids_all <- omic_features[[id_col]]

  sample_z <- scale(t(numeric_df))
  sample_outliers <- abs(sample_z) > 3
  sample_outliers[is.na(sample_outliers)] <- FALSE
  colnames(sample_outliers) <- sample_ids_all

  n_outliers_sample <- colSums(sample_outliers, na.rm = TRUE)
  max_abs_z_sample  <- suppressWarnings(apply(abs(sample_z), 2, max, na.rm = TRUE))

  sample_severity <- data.frame(
    sample     = colnames(sample_outliers),
    n_outliers = n_outliers_sample,
    max_abs_z  = max_abs_z_sample
  )
  sample_severity <- sample_severity[sample_severity$n_outliers > 0, ]
  sample_severity <- sample_severity[order(-sample_severity$max_abs_z), ]

  sample_medians <- sapply(seq_len(nrow(numeric_df)), function(i) {
    median(as.numeric(numeric_df[i, ]), na.rm = TRUE)
  })
  names(sample_medians) <- sample_ids_all

  outlier_sample_ids <- names(n_outliers_sample[n_outliers_sample >= 1])

  # order samples from min-max median intensity
  outlier_sample_ids <- outlier_sample_ids[order(sample_medians[outlier_sample_ids])]

  cat("Total samples with outliers (>3 SD):", length(outlier_sample_ids), "\n")

  # skip forest plots if no outliers detected
  if (length(outlier_sample_ids) == 0) {
    cat("No samples with outliers detected. Forest plots skipped.\n\n")
    all_batch_stats <- data.frame()
  } else {

  # batch the samples
    batch_size <- if (length(outlier_sample_ids) <= 100) length(outlier_sample_ids) else 60

    sample_ids  <- outlier_sample_ids
    n_samples   <- length(sample_ids)
    n_batches   <- ceiling(n_samples / batch_size)
    all_batch_stats <- list()

    cat("\nGenerating forest plots...\n")
    outlier_progress <- txtProgressBar(min = 0, max = n_batches, style = 3)

    for (b in seq_len(n_batches)) {

      start_idx <- (b - 1) * batch_size + 1
      end_idx   <- min(b * batch_size, n_samples)
      batch_ids <- sample_ids[start_idx:end_idx]

      batch_stats <- lapply(batch_ids, function(samp) {
        row_idx <- which(sample_ids_all == samp)
        values  <- as.numeric(numeric_df[row_idx, ])
        values  <- values[!is.na(values)]

        outlier_flag       <- sample_outliers[, samp]
        non_outlier_values <- values[!outlier_flag]
        if (length(non_outlier_values) == 0) non_outlier_values <- values

      c(
        median = median(values),
        q1     = quantile(values, 0.25),
        q3     = quantile(values, 0.75),
        wmin   = min(non_outlier_values),
        wmax   = max(non_outlier_values)
      )
    })

    batch_df <- data.frame(
      sample = batch_ids,
      median = sapply(batch_stats, `[`, "median"),
      q1     = sapply(batch_stats, `[`, "q1.25%"),
      q3     = sapply(batch_stats, `[`, "q3.75%"),
      wmin   = sapply(batch_stats, `[`, "wmin"),
      wmax   = sapply(batch_stats, `[`, "wmax")
    )

      all_batch_stats[[b]] <- batch_df

      batch_full_range <- range(
        sapply(batch_ids, function(samp) {
          row_idx <- which(sample_ids_all == samp)
          as.numeric(numeric_df[row_idx, ])
        }),
        na.rm = TRUE
      )

      par(mfrow = c(1, 1), mar = c(8, 4, 3, 2))
      x_pos      <- seq_len(nrow(batch_df))
      box_width  <- 0.1

      plot(x_pos, batch_df$median,
          ylim = batch_full_range,
          xlim = c(0.5, nrow(batch_df) + 0.5),
          xaxt = "n", xlab = "", ylab = "Relative Intensity",
          type = "n",
          main = paste0("Sample Intensity with Outliers Summary (Plot ",
                       b, " of ", n_batches, ")"))

    # whiskers
      segments(x_pos, batch_df$wmin, x_pos, batch_df$q1, col = "black", lwd = 1)
      segments(x_pos, batch_df$q3,   x_pos, batch_df$wmax, col = "black", lwd = 1)

    # IQR box
      rect(
        xleft   = x_pos - box_width,
        ybottom = batch_df$q1,
        xright  = x_pos + box_width,
        ytop    = batch_df$q3,
        col     = "white",
        border  = "black"
      )

    # median point on top of the box
      points(x_pos, batch_df$median, pch = 16, col = "black", cex = 0.8)

    # x-axis labels
      axis(1, at = x_pos, labels = batch_df$sample, las = 2, cex.axis = 0.6)

    # overlay outliers as open black circles
      for (i in seq_along(batch_ids)) {
        samp         <- batch_ids[i]
        row_idx      <- which(sample_ids_all == samp)
        values       <- as.numeric(numeric_df[row_idx, ])
        outlier_flag <- sample_outliers[, samp]
        if (any(outlier_flag, na.rm = TRUE)) {
          points(rep(i, sum(outlier_flag, na.rm = TRUE)), values[outlier_flag],
                col = "black", pch = 1, cex = 1.2)
        }
      }

      setTxtProgressBar(outlier_progress, b)
    }

    close(outlier_progress)
    par(mfrow = c(1, 1))
    all_batch_stats <- do.call(rbind, all_batch_stats)
  }

  # Step 2: PCA - Variance Explained, Scree, Loadings, Batch Assessment

  # drop incomplete rows
  complete_rows <- complete.cases(numeric_df)
  pca_input <- numeric_df[complete_rows, ]

  cat("Samples used in PCA (complete cases):", nrow(pca_input),
      "out of", nrow(numeric_df), "\n")

  # check variance on those rows
  col_variances <- apply(pca_input, 2, var, na.rm = TRUE)
  pca_input <- pca_input[, col_variances > 0 & !is.na(col_variances), drop = FALSE]

  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)

  # variance explained table for top 5 PCs
  var_summary <- summary(pca_result)$importance
  n_pcs_available <- ncol(pca_result$rotation)

  # show up to 5 or less if necessary
  pca_variance <- data.frame(
    PC                  = colnames(var_summary),
    proportion_variance = round(var_summary["Proportion of Variance", ], 4),
    cumulative_variance = round(var_summary["Cumulative Proportion", ], 4)
  )
  n_pcs_variance <- min(5, n_pcs_available)
  pca_variance   <- head(pca_variance, n_pcs_variance)

  cat("\nVariance Explained by Top", n_pcs_variance, "Principal Components:\n")
  print(pca_variance)

  # up to top 10 absolute loadings per PC, for up to top 5 PCs
  n_pcs_show <- min(5, n_pcs_available)
  top_n_loadings <- 10
  loadings <- pca_result$rotation

  loading_table <- do.call(rbind, lapply(seq_len(n_pcs_show), function(i) {
    pc_loadings <- loadings[, i]
    n_top       <- min(top_n_loadings, nrow(loadings))
    top_idx <- order(abs(pc_loadings), decreasing = TRUE)[seq_len(n_top)]
    data.frame(
      PC      = paste0("PC", i),
      feature = rownames(loadings)[top_idx],
      loading = round(pc_loadings[top_idx], 4),
      row.names = NULL
    )
  }))
  rownames(loading_table) <- NULL

  cat("\nTop", top_n_loadings, "Absolute Loadings for Top", n_pcs_show, "PCs:\n")
  cat("Showing first 5 rows. Full table available in returned object.\n")
  print(head(loading_table, 5))

  n_pcs_label <- min(2, n_pcs_available)
  var_explained <- round(var_summary[2, seq_len(n_pcs_label)] * 100, 1)

  # batch setup and run ANOVA if batch is defined by user
  anova_pc1            <- NULL
  anova_pc2            <- NULL
  p_pc1                <- NULL
  p_pc2                <- NULL
  point_colors         <- "steelblue"
  batch_levels         <- NULL
  batch_colors_palette <- NULL

  if (!is.null(batch)) {
    batch_pca    <- as.factor(batch[complete_rows])
    batch_levels <- levels(batch_pca)

    # Okabe-Ito palette
    okabe_ito <- c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#999999"
    )

    # fallback to hcl.colors if > 8 batches
    if (length(batch_levels) <= 8) {
      batch_colors_palette <- okabe_ito[seq_along(batch_levels)]
    } else {
      batch_colors_palette <- hcl.colors(length(batch_levels), palette = "Dark 3")
    }

    point_colors <- batch_colors_palette[match(batch_pca, batch_levels)]

    # build dataframe for ANOVA
    pca_df_batch <- data.frame(
      PC1   = pca_result$x[, 1],
      PC2   = pca_result$x[, 2],
      batch = batch_pca
    )

    # ANOVA: PC1 ~ batch
    anova_pc1 <- aov(PC1 ~ batch, data = pca_df_batch)
    p_pc1     <- summary(anova_pc1)[[1]][["Pr(>F)"]][1]

    # ANOVA: PC2 ~ batch
    anova_pc2 <- aov(PC2 ~ batch, data = pca_df_batch)
    p_pc2     <- summary(anova_pc2)[[1]][["Pr(>F)"]][1]

    cat("\nBatch Effect Assessment (ANOVA):\n")
    if (!is.na(p_pc1)) {
      if (p_pc1 < 0.05) {
        cat("SIGNIFICANT batch effect on PC1 (p =", signif(p_pc1, 3), ")\n")
      } else {
        cat("No significant batch effect on PC1 (p =", signif(p_pc1, 3), ")\n")
      }
    }
    if (!is.na(p_pc2)) {
      if (p_pc2 < 0.05) {
        cat("SIGNIFICANT batch effect on PC2 (p =", signif(p_pc2, 3), ")\n")
      } else {
        cat("No significant batch effect on PC2 (p =", signif(p_pc2, 3), ")\n")
      }
    }
  }

  if (!is.null(batch)) {

  # stacked boxplots of PC1 and PC2 by batch
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

  # PC1 boxplot
    boxplot(PC1 ~ batch, data = pca_df_batch,
          main  = paste0("PC1 by Batch (p = ", signif(p_pc1, 3), ")"),
          xlab  = "Batch",
          ylab  = "PC1",
          col   = batch_colors_palette,
          border = "black")

    stripchart(PC1 ~ batch, data = pca_df_batch,
             vertical = TRUE,
             method   = "jitter",
             add      = TRUE,
             pch      = 16,
             col      = "black",
             cex      = 0.6)

  # PC2 boxplot
    boxplot(PC2 ~ batch, data = pca_df_batch,
          main  = paste0("PC2 by Batch (p = ", signif(p_pc2, 3), ")"),
          xlab  = "Batch",
          ylab  = "PC2",
          col   = batch_colors_palette,
          border = "black")

    stripchart(PC2 ~ batch, data = pca_df_batch,
             vertical = TRUE,
             method   = "jitter",
             add      = TRUE,
             pch      = 16,
             col      = "black",
             cex      = 0.6)
  }
  # scree plot + PCA scatterplot
    par(mfrow = c(1, 2))

    n_pcs_scree <- min(15, n_pcs_available)

    barplot(var_summary["Proportion of Variance", seq_len(n_pcs_scree)] * 100,
          names.arg = colnames(var_summary)[seq_len(n_pcs_scree)],
          las  = 2,
          col  = "steelblue",
          main = "Scree Plot",
          ylab = "% Variance Explained")

    plot(pca_result$x[, 1], pca_result$x[, 2],
       col  = point_colors,
       pch  = 16,
       xlab = paste0("PC1 (", var_explained[1], "%)"),
       ylab = paste0("PC2 (", var_explained[2], "%)"),
       main = if (!is.null(batch)) "PCA (colored by batch)" else "PC1 vs PC2")

    # Step 3: Correlation Blocks via Hierarchical Clustering
    cat("\nComputing correlation matrix...\n")

    X <- as.matrix(numeric_df)

    # remove zero-variance columns before correlation
    col_vars     <- apply(X, 2, var, na.rm = TRUE)
    zero_var_cols <- sum(col_vars == 0 | is.na(col_vars))
    X <- X[, col_vars > 0 & !is.na(col_vars), drop = FALSE]

    if (zero_var_cols > 0) {
      cat(zero_var_cols, "column(s) with zero standard deviation removed",
          "before computing correlation matrix.\n")
    }

    n_feat     <- ncol(X)
    chunk_size <- 50
    n_chunks   <- ceiling(n_feat / chunk_size)

    cor_mat <- matrix(NA, n_feat, n_feat,
                      dimnames = list(colnames(X), colnames(X)))

    cor_progress <- txtProgressBar(min = 0, max = n_chunks, style = 3)

    for (c in seq_len(n_chunks)) {
      start_idx <- (c - 1) * chunk_size + 1
      end_idx   <- min(c * chunk_size, n_feat)
      col_idx   <- start_idx:end_idx

      cor_mat[, col_idx] <- cor(
        X, X[, col_idx, drop = FALSE],
        use = "pairwise.complete.obs", method = "spearman"
      )

      setTxtProgressBar(cor_progress, c)
    }

    close(cor_progress)

    # check for any remaining NAs in the corr matrix
    na_cols <- apply(cor_mat, 1, function(x) any(is.na(x)))
    n_na_cols <- sum(na_cols)

    if (n_na_cols > 0) {
      cat(n_na_cols, "column(s) produced NA correlations",
          "and were removed before clustering.\n")
      cor_mat <- cor_mat[!na_cols, !na_cols]
    }

    # convert correlation to distance, then cluster
    d  <- as.dist(1 - cor_mat)
    hc <- hclust(d, method = "complete")

    blocks <- cutree(hc, h = 1 - rho_thresh)
    block_df <- data.frame(
      feature = names(blocks),
      block   = as.integer(blocks)
    )

    block_sizes <- aggregate(feature ~ block, data = block_df, FUN = length)
    names(block_sizes) <- c("block", "n_features")
    block_sizes <- block_sizes[order(-block_sizes$n_features), ]

    keep_blocks <- block_sizes$block[block_sizes$n_features >= min_block_size]
    block_df$block_kept <- block_df$block %in% keep_blocks

    kept_blocks <- block_df[block_df$block_kept, ]
    kept_block_features <- lapply(split(kept_blocks$feature,
                                        kept_blocks$block), identity)
    kept_block_summary <- data.frame(
      block      = names(kept_block_features),
      n_features = as.integer(sapply(kept_block_features, length)),
      row.names  = NULL
    )
    kept_block_summary <- kept_block_summary[order(-kept_block_summary$n_features), ]
    rownames(kept_block_summary) <- NULL

    if (nrow(kept_block_summary) == 0) {
      cat("\nNo correlation blocks detected meeting the specified criteria ",
          "(rho >=", rho_thresh, ", min size", min_block_size, ").\n")
    } else {
      cat("\nCorrelation Blocks (rho >=", rho_thresh, ", min size",
          min_block_size, "):\n")
      cat("Showing first 5 of", nrow(kept_block_summary),
          "blocks. Full table available in returned object.\n")
      print(head(kept_block_summary, 5))
    }

  invisible(list(
    sample_severity     = sample_severity,
    sample_outliers     = sample_outliers,
    batch_stats         = all_batch_stats,
    pca_result          = pca_result,
    pca_variance        = pca_variance,
    pca_loadings        = loading_table,
    anova_pc1           = anova_pc1,
    anova_pc2           = anova_pc2,
    p_value_pc1         = p_pc1,
    p_value_pc2         = p_pc2,
    correlation_blocks  = kept_block_features,
    correlation_summary = kept_block_summary
  ))
  
  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
  
}
