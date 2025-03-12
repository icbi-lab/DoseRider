#' Perform DoseRider Analysis Using Linear Mixed Models (LMMs) or Generalized Additive Mixed Models (GAMMs)
#'
#' This function performs DoseRider analysis on gene expression data, applying either Linear Mixed Models (LMMs) or
#' Generalized Additive Mixed Models (GAMMs) to each gene set defined in the gene matrix transposed (GMT) format.
#' It evaluates dose-response relationships in the context of gene sets and calculates various model metrics,
#' significance, Benchmark Dose (BMD), and smoothing predictions.
#'
#' @param se A `SummarizedExperiment` object or a matrix/data frame containing gene expression data.
#' @param gmt A list of gene sets, where each entry contains a list of gene names.
#' @param dose_col A string specifying the column name representing dose information.
#' @param sample_col A string specifying the column name representing sample information.
#' @param covariates Optional. A character vector specifying additional covariate column(s) in `se`. Default is an empty vector.
#' @param omic A string specifying the type of omics data. Default is `"rnaseq"`.
#' @param minGSsize An integer specifying the minimum gene set size for analysis. Default is `5`.
#' @param maxGSsize An integer specifying the maximum gene set size for analysis. Default is `300`.
#' @param method A string specifying the method for multiple testing adjustment. Default is `"fdr"`.
#' @param modelType A string specifying the type of model to fit for each gene set. Options:
#'   \itemize{
#'     \item `"LMM"`: Linear Mixed Model
#'     \item `"GAMM"`: Generalized Additive Mixed Model
#'   }
#'   Default is `"LMM"`.
#' @param FilterPathway Logical. Whether to apply pathway-level filtering. Default is `FALSE`.
#' @param filters A character vector specifying the filtering methods to apply. Options:
#'   \itemize{
#'     \item `"pca"`: Principal Component Analysis (PCA)-based filtering to detect antagonistic patterns.
#'     \item `"variance"`: Variance-based filtering to remove gene sets with low expression variability.
#'     \item `"correlation"`: Correlation-based filtering to remove pathways with highly antagonistic gene expression patterns.
#'   }
#' @param filter_params A named list of filter-specific parameters. Available options:
#'   \itemize{
#'     \item **For PCA Filtering (`"pca"`)**:
#'       \itemize{
#'         \item `pca_threshold` (numeric): Variance threshold for PC1 to filter pathways. Default: `0.6`.
#'         \item `loading_threshold` (numeric): Threshold to detect antagonistic patterns based on PCA loadings. Default: `0.6`.
#'         \item `antagonist_threshold` (numeric): Defines when a pathway is considered antagonistic. Default: `2.0`.
#'       }
#'     \item **For Variance Filtering (`"variance"`)**:
#'       \itemize{
#'         \item `variance_percentile` (numeric): Percentile threshold for gene variance filtering. Default: `0.2` (removes genes in the lowest 20% of variance).
#'       }
#'     \item **For Correlation Filtering (`"correlation"`)**:
#'       \itemize{
#'         \item `correlation_threshold` (numeric): Proportion of negatively correlated genes required to exclude a pathway. Default: `0.5`.
#'         \item `negative_corr_cutoff` (numeric): Correlation coefficient threshold below which genes are considered antagonistic. Default: `-0.3`.
#'       }
#'   }
#' @param log_transform Logical. Whether to log10 transform the dose values. Default is `FALSE`.
#' @param spline_knots An integer specifying the number of internal knots for the spline function. Default is `2`.
#'                     The number of knots should be lower than the number of unique dose values.
#' @param knot_method A string specifying the method for selecting spline knots. Options:
#'   \itemize{
#'     \item `"quantile"` (default): Selects knots based on percentiles of the dose distribution.
#'     \item `"geometric"`: Places knots at logarithmically spaced intervals.
#'     \item `"manual"`: Uses user-defined knots from the `manual_knots` parameter.
#'   }
#' @param manual_knots A numeric vector specifying custom knot positions (used only if `knot_method = "manual"`).
#'                     If `NULL`, this option is ignored.
#'
#' @return A list containing results for each gene set, including:
#'   \itemize{
#'     \item Model metrics (e.g., AIC, BIC)
#'     \item Significance tests (p-values and adjusted p-values)
#'     \item Benchmark Dose (BMD) calculations
#'     \item Smoothed dose-response predictions
#'   }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("SummarizedExperiment")
#'
#' # Define a small gene set
#' gmt <- list(geneSet1 = list(genes = c("gene1", "gene2")))
#'
#' # Run DoseRider with default settings (LMM model)
#' results <- DoseRider(se, gmt, dose_col = "dose", sample_col = "sample",
#'                      omic = "rnaseq", modelType = "LMM")
#'
#' # Run DoseRider using a GAMM model with geometric spline knots
#' results_gamm <- DoseRider(se, gmt, dose_col = "dose", sample_col = "sample",
#'                           modelType = "GAMM", knot_method = "geometric", spline_knots = 3)
#'
#' # Run DoseRider with manual knot selection
#' manual_knot_positions <- c(0.01, 0.1, 1, 10)
#' results_manual <- DoseRider(se, gmt, dose_col = "dose", sample_col = "sample",
#'                             modelType = "LMM", knot_method = "manual", manual_knots = manual_knot_positions)
#' }
#'
#' @importFrom stats p.adjust
#' @import utils
#' @export
process_gene_set <- function(se, dose_col, sample_col, omic, gmt, i, minGSsize = 5,
                             maxGSsize = 300, covariates = c(), modelType = "LMM",
                             FilterPathway = FALSE,
                             filters = c("correlation"),
                             filter_params = list(),
                             models = c("linear", "non_linear_fixed","non_linear_mixed"),
                             log_transform = F,
                             spline_knots = 2, knot_method = "quantile",
                             manual_knots = NULL) {
  # Helper function to summarize a model
  is_fitted_model <- function(model) {
    if(inherits(model, c("lmerMod", "glmerMod"))) {
      return(TRUE)
    } else if(inherits(model, c("gam", "bam"))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }

  geneset <- gmt[[i]]$genes
  # Prepare data for the gene set   return(list(long_df = long_df, spline_info = spline_info))
  data_info <- suppressWarnings(prepare_data(se=se, geneset=geneset,
                                dose_col=dose_col, sample_col=sample_col,
                                omic=omic, log_transform=log_transform,
                                spline_knots = spline_knots, knot_method = knot_method,
                                manual_knots = manual_knots))
  long_df <- data_info$long_df
  spline_info <- data_info$spline_info
  if (log_transform) {
    dose_col <- paste0("log_", dose_col)

  }
  # Check gene set size criteria
  if (is.null(long_df) || length(unique(long_df$gene)) < minGSsize || length(unique(long_df$gene)) > maxGSsize) {
    return(NULL)
  }

  # Apply filtering if enabled
  if (FilterPathway) {
    should_keep_pathway <- apply_pathway_filters(long_df, filters = filters, filter_params = filter_params)
    if (!should_keep_pathway) {
      return(NULL)
    }
  }

  # Fit models according to user-specified models and always include the base (null) model
  fitted_models <- list(null = NULL, linear = NULL, non_linear_fixed = NULL, non_linear_mixed = NULL)

  if (modelType == "LMM") {
    # Always fit the base (null) model
    null_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "null", omic)
    null_results <- suppressWarnings(fit_lmm(null_formula, long_df, omic))
    fitted_models$null <- null_results

    # Fit the models specified by the user
    for (model in models) {
      if (model == "linear") {
        linear_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "linear", omic,spline_knots)
        fitted_models$linear <- suppressWarnings(fit_lmm(linear_formula, long_df, omic))
      } else if (model == "non_linear_fixed") {
        non_linear_fixed_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "non_linear_fixed", omic,spline_knots)
        fitted_models$non_linear_fixed <- suppressWarnings(fit_lmm(non_linear_fixed_formula, long_df, omic))
      } else if (model == "non_linear_mixed") {
        non_linear_mixed_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "non_linear_mixed", omic,spline_knots)
        fitted_models$non_linear_mixed <- suppressWarnings(fit_lmm(non_linear_mixed_formula, long_df, omic))
      }
    }
  } else if (modelType == "GAMM") {
    null_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "null", omic)
    fitted_models$null <- suppressWarnings(fit_gam(null_formula, long_df, omic))

    for (model in models) {
      if (model == "linear") {
        linear_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "linear", omic)
        fitted_models$linear <- suppressWarnings(fit_gam(linear_formula, long_df, omic))
      } else if (model == "non_linear_fixed") {
        non_linear_fixed_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "non_linear_fixed", omic)
        fitted_models$non_linear_fixed <- suppressWarnings(fit_gam(non_linear_fixed_formula, long_df, omic))
      } else if (model == "non_linear_mixed") {
        non_linear_mixed_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "non_linear_mixed", omic)
        fitted_models$non_linear_mixed <- suppressWarnings(fit_gam(non_linear_mixed_formula, long_df, omic))
      }
    }
  }
  is_fitted <- sapply(fitted_models, is_fitted_model)
  filter_fitted_models <- fitted_models[is_fitted]

  if (length(filter_fitted_models) > 1) {
    p_value_list <- compare_all_models(fitted_models$nul,
                                       fitted_models$linear,
                                       fitted_models$non_linear_fixed,
                                       fitted_models$non_linear_mixed,
                                       modelType)

    best_model_AICc <- select_best_model(fitted_models, p_value_list)
    best_model <- fitted_models[[best_model_AICc]]
  } else {
    p_value_list <- list("p_value_linear"= NA, "p_value_non_linear"= NA)
    best_model_AICc <- names(fitted_models)
    best_model <- fitted_models[[1]]
  }

  # Compute metrics for null and cubic models
  null_metrics <- if (modelType=="LMM") compute_metrics_lmm(fitted_models$null) else compute_metrics_gamm(fitted_models$null)
  linear_metrics <- if (modelType=="LMM") compute_metrics_lmm(fitted_models$linear) else compute_metrics_gamm(fitted_models$linear)
  non_linear_fixed_metrics <- if (modelType=="LMM") compute_metrics_lmm(fitted_models$non_linear_fixed) else compute_metrics_gamm(fitted_models$non_linear_fixed)
  non_linear_mixed_metrics <- if (modelType=="LMM") compute_metrics_lmm(fitted_models$non_linear_mixed) else compute_metrics_gamm(fitted_models$non_linear_mixed)

  # Calculate smoothing values and BMD for the best model
  if (best_model_AICc != "null") {
    # Compute the smooth values for the best model
    smooth_values <- smooth_pathway_trend(best_model, long_df, dose_col, sample_col, omic, TRUE, covariates, dose_points = 50, spline_info = spline_info)

    # Extract random effects based on the model type
    random_effect <- if (modelType == "LMM") extract_random_effects_lmm(best_model, dose_col) else extract_random_effects_gamm(best_model, dose_col)

    # Add predictions to long_df
    long_df$predictions <- predict(best_model, newdata = long_df)

    # Initialize results list
    cluster_specific_results <- list()

    # Step 1: Compute the overall results for all genes (AllGenes)
    smooth_cluster_all <- smooth_values
    derivative_cluster_all <- compute_derivatives(smooth_cluster_all, dose_col)
    bmd_cluster_all <- compute_bmd_from_main_trend(smooth_cluster_all, dose_col, z = 1)

    # Store results for AllGenes
    cluster_specific_results[["AllGenes"]] <- list(
      Derivative = derivative_cluster_all,
      BMD = bmd_cluster_all
    )

    # Assign "AllGenes" to all genes in the smooth_values dataframe
    all_genes <- unique(smooth_values$gene)
    cluster <- rep("AllGenes", length(all_genes))
    names(cluster) <- all_genes

    # Step 2: Compute the clustering results
    # Determine the number of unique genes
    num_genes <- length(unique(long_df$gene))

    # Set max clusters to either the number of unique genes - 1 or 5, whichever is smaller
    max_clusters <- min(5, num_genes - 1)

    # Compute the optimal number of clusters based on silhouette analysis
    optimal_clusters_silhouette <- optimal_clusters_silhouette(smooth_values, dose_col, max_clusters = max_clusters)

    # Extract cluster assignments and number of optimal clusters
    cluster <- optimal_clusters_silhouette$Cluster
    n_cluster <- optimal_clusters_silhouette$OptimalClusters

    # Compute mean trend, BMD, and derivatives for each cluster
    for (j in seq_len(n_cluster)) {
      cluster_genes <- names(cluster[cluster == j])
      smooth_cluster <- smooth_values[smooth_values$gene %in% cluster_genes, ]

      # Compute the derivatives and BMD for each cluster
      derivative_cluster <- compute_derivatives(smooth_cluster, dose_col)
      bmd_cluster <- compute_bmd_from_main_trend(smooth_cluster, dose_col, z = 1)

      # Store cluster-specific results
      cluster_specific_results[[paste("Cluster", j)]] <- list(
        Derivative = derivative_cluster,
        BMD = bmd_cluster
      )
    }

  } else {
    # If there is no valid model, initialize empty results
    smooth_values <- NA
    random_effect <- NA
    cluster <- NA
    n_cluster <- NA
    cluster_specific_results <- NA
  }


  # Compile results into a list
  geneset_results <- list(
    Geneset = gmt[[i]]$pathway,
    Geneset_Size = length(geneset),
    Genes = length(unique(long_df$gene)),
    null_AIC = null_metrics$AIC,
    null_AICc = null_metrics$AICc,
    null_BIC = null_metrics$BIC,
    null_df = null_metrics$edf,
    null_df = null_metrics$ICC,
    linear_AIC = linear_metrics$AIC,
    linear_AICc = linear_metrics$AICc,
    linear_BIC = linear_metrics$BIC,
    linear_df = linear_metrics$edf,
    linear_ICC = linear_metrics$ICC,
    non_linear_fixed_AIC = non_linear_fixed_metrics$AIC,
    non_linear_fixed_AICc = non_linear_fixed_metrics$AICc,
    non_linear_fixed_BIC = non_linear_fixed_metrics$BIC,
    non_linear_fixed_df = non_linear_fixed_metrics$edf,
    non_linear_fixed_ICC = non_linear_fixed_metrics$ICC,
    non_linear_mixed_AIC = non_linear_mixed_metrics$AIC,
    non_linear_mixed_AICc = non_linear_mixed_metrics$AICc,
    non_linear_mixed_BIC = non_linear_mixed_metrics$BIC,
    non_linear_mixed_df = non_linear_mixed_metrics$edf,
    non_linear_mixed_ICC = non_linear_mixed_metrics$ICC,
    p_value_linear = p_value_list$p_value_linear,
    p_value_non_linear_fixed = p_value_list$p_value_non_linear_fixed,
    p_value_non_linear_mixed = p_value_list$p_value_non_linear_mixed,

    best_model = best_model_AICc,
    best_model_pvalue = ifelse(best_model_AICc == "null", NA,
                               ifelse(best_model_AICc == "linear", p_value_list$p_value_linear,
                                      ifelse(best_model_AICc == "non_linear_fixed", p_value_list$p_value_non_linear_fixed,
                                             ifelse(best_model_AICc == "non_linear_mixed", p_value_list$p_value_non_linear_mixed, NA)))),


    best_model_aicc = ifelse(best_model_AICc == "null", null_metrics$AICc,
                             ifelse(best_model_AICc == "linear", linear_metrics$AICc,
                                    ifelse(best_model_AICc == "non_linear_fixed", non_linear_fixed_metrics$AICc,
                                           ifelse(best_model_AICc == "non_linear_mixed", non_linear_mixed_metrics$AICc, NA)))),

    best_model_bic = ifelse(best_model_AICc == "null", null_metrics$BIC,
                            ifelse(best_model_AICc == "linear", linear_metrics$BIC,
                                   ifelse(best_model_AICc == "non_linear_fixed", non_linear_fixed_metrics$BIC,
                                          ifelse(best_model_AICc == "non_linear_mixed", non_linear_mixed_metrics$BIC, NA)))),
    Smooth_Predictions = list(smooth_values),
    #Smooth_Predictions_Pathway = list(smooth_pathway),
    Raw_Values = list(long_df),
    #BMD = bmd,
    #TCD = derivative,
    random_effect = random_effect,
    ClusterAssignments = cluster,
    OptimalClusters = n_cluster,
    ClusterSpecificResults = cluster_specific_results,
    dose_col = dose_col,
    log_transform = log_transform
  )

  return(geneset_results)
}

#' Perform DoseRider Analysis Using Linear Mixed Models (LMMs) or Generalized Additive Mixed Models (GAMMs)
#'
#' This function performs DoseRider analysis on gene expression data, applying either Linear Mixed Models (LMMs)
#' or Generalized Additive Mixed Models (GAMMs) to each gene set defined in the gene matrix transposed (GMT) format.
#' It evaluates dose-response relationships in the context of gene sets and calculates various model metrics,
#' significance, and smoothing predictions.
#'
#' @param se A `SummarizedExperiment` object or a matrix/data frame containing gene expression data.
#' @param gmt A list of gene sets, where each entry contains a list of gene names.
#' @param dose_col A string specifying the column name representing dose information.
#' @param sample_col A string specifying the column name representing sample information.
#' @param covariates Optional. A character vector specifying additional covariate column(s) in `se`. Default is an empty vector.
#' @param omic A string specifying the type of omics data. Default is `"rnaseq"`.
#' @param minGSsize An integer specifying the minimum gene set size for analysis. Default is `5`.
#' @param maxGSsize An integer specifying the maximum gene set size for analysis. Default is `300`.
#' @param method A string specifying the method for multiple testing adjustment. Default is `"fdr"`.
#' @param modelType A string specifying the type of model to fit for each gene set. Options:
#'   \itemize{
#'     \item `"LMM"`: Linear Mixed Model
#'     \item `"GAMM"`: Generalized Additive Mixed Model
#'   }
#'   Default is `"LMM"`.
#' @param FilterPathway Logical. Whether to apply pathway-level filtering. Default is `FALSE`.
#' @param filters A character vector specifying the filtering methods to apply. Options:
#'   \itemize{
#'     \item `"pca"`: Principal Component Analysis (PCA)-based filtering to detect antagonistic patterns.
#'     \item `"variance"`: Variance-based filtering to remove gene sets with low expression variability.
#'     \item `"correlation"`: Correlation-based filtering to remove pathways with highly antagonistic gene expression patterns.
#'   }
#' @param filter_params A named list of filter-specific parameters. Available options:
#'   \itemize{
#'     \item **For PCA Filtering (`"pca"`)**:
#'       \itemize{
#'         \item `pca_threshold` (numeric): Variance threshold for PC1 to filter pathways. Default: `0.6`.
#'         \item `loading_threshold` (numeric): Threshold to detect antagonistic patterns based on PCA loadings. Default: `0.6`.
#'         \item `antagonist_threshold` (numeric): Defines when a pathway is considered antagonistic. Default: `2.0`.
#'       }
#'     \item **For Variance Filtering (`"variance"`)**:
#'       \itemize{
#'         \item `variance_percentile` (numeric): Percentile threshold for gene variance filtering. Default: `0.2` (removes genes in the lowest 20% of variance).
#'       }
#'     \item **For Correlation Filtering (`"correlation"`)**:
#'       \itemize{
#'         \item `correlation_threshold` (numeric): Proportion of negatively correlated genes required to exclude a pathway. Default: `0.5`.
#'         \item `negative_corr_cutoff` (numeric): Correlation coefficient threshold below which genes are considered antagonistic. Default: `-0.3`.
#'       }
#'   }
#' @param log_transform Logical. Whether to log10 transform the dose values. Default is `FALSE`.
#' @param spline_knots An integer specifying the number of internal knots for the spline function. Default is `2`.
#'                     The number of knots should be lower than the number of unique dose values.
#' @param knot_method A string specifying the method for selecting spline knots. Options:
#'   \itemize{
#'     \item `"quantile"` (default): Selects knots based on percentiles of the dose distribution.
#'     \item `"geometric"`: Places knots at logarithmically spaced intervals.
#'     \item `"manual"`: Uses user-defined knots from the `manual_knots` parameter.
#'   }
#' @param manual_knots A numeric vector specifying custom knot positions (used only if `knot_method = "manual"`).
#'                     If `NULL`, this option is ignored.
#'
#'
#' @return A list containing results for each gene set including various metrics, p-values,
#' and adjusted p-values.
#'
#' @examples
#' \dontrun{
#' data("SummarizedExperiment")
#' gmt <- list(geneSet1 = list(genes = c("gene1", "gene2")))
#' results <- DoseRider(se, gmt, "dose", "sample", "covariate", "rnaseq", modelType = "GAMM")
#' }
#'
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar
#' @import utils
#' @import progress
#' @export
DoseRider <- function(se, gmt, dose_col = "dose", sample_col = "sample",
                      covariates = c(), omic = "rnaseq", minGSsize = 5,
                      maxGSsize = 300, method = "fdr", modelType = "LMM",
                      FilterPathway = FALSE,
                      filters = c("correlation"),
                      filter_params = list(),
                      log_transform = F,
                      models = c("linear", "non_linear_fixed","non_linear_mixed"),
                      spline_knots = 2, knot_method = "quantile",
                      manual_knots = NULL) {

  # Validate input data
  validate_input_doserider(se, dose_col, sample_col, covariates)

  # Initialize results list
  results <- list()
  total_gene_sets <- length(gmt)
  pb <- txtProgressBar(min = 0, max = total_gene_sets, style = 3)

  # Process each gene set
  for (i in seq_along(gmt)) {
    setTxtProgressBar(pb, i)
    geneset_results <- suppressMessages(suppressWarnings(process_gene_set(
      se = se, dose_col = dose_col, sample_col = sample_col, omic = omic, gmt = gmt,
      i = i, minGSsize = minGSsize, maxGSsize = maxGSsize, covariates = covariates,
      modelType = modelType, FilterPathway =FilterPathway,
      filters = filters,
      filter_params = filter_params,
      log_transform = log_transform,
      models = models
    )))
    if (!is.null(geneset_results)) {
      results[[gmt[[i]]$pathway]] <- geneset_results
    }
  }

  close(pb)

  # Adjust p-values for multiple testing
  results <- adjust_pvalues_doserider_result(results = results, method = method)

  class(results) <- "DoseRider"
  return(results)
}


#' Perform DoseRider analysis in parallel using multiple cores.
#'
#' This function conducts DoseRider analysis in parallel, applying either Linear Mixed Models (LMMs) or Generalized
#' Additive Mixed Models (GAMMs) to each gene set in the gene matrix transposed (GMT) format. It evaluates dose-response
#' relationships in gene sets, computing various model metrics and significance.
#'
#' @param se A `SummarizedExperiment` object or a matrix/data frame containing gene expression data.
#' @param gmt A list of gene sets, where each entry contains a list of gene names.
#' @param dose_col A string specifying the column name representing dose information.
#' @param sample_col A string specifying the column name representing sample information.
#' @param covariates Optional. A character vector specifying additional covariate column(s) in `se`. Default is an empty vector.
#' @param omic A string specifying the type of omics data. Default is `"rnaseq"`.
#' @param minGSsize An integer specifying the minimum gene set size for analysis. Default is `5`.
#' @param maxGSsize An integer specifying the maximum gene set size for analysis. Default is `300`.
#' @param method A string specifying the method for multiple testing adjustment. Default is `"fdr"`.
#' @param modelType A string specifying the type of model to fit for each gene set. Options:
#'   \itemize{
#'     \item `"LMM"`: Linear Mixed Model
#'     \item `"GAMM"`: Generalized Additive Mixed Model
#'   }
#'   Default is `"LMM"`.
#' @param FilterPathway Logical. Whether to apply pathway-level filtering. Default is `FALSE`.
#' @param filters A character vector specifying the filtering methods to apply. Options:
#'   \itemize{
#'     \item `"pca"`: Principal Component Analysis (PCA)-based filtering to detect antagonistic patterns.
#'     \item `"variance"`: Variance-based filtering to remove gene sets with low expression variability.
#'     \item `"correlation"`: Correlation-based filtering to remove pathways with highly antagonistic gene expression patterns.
#'   }
#' @param filter_params A named list of filter-specific parameters. Available options:
#'   \itemize{
#'     \item **For PCA Filtering (`"pca"`)**:
#'       \itemize{
#'         \item `pca_threshold` (numeric): Variance threshold for PC1 to filter pathways. Default: `0.6`.
#'         \item `loading_threshold` (numeric): Threshold to detect antagonistic patterns based on PCA loadings. Default: `0.6`.
#'         \item `antagonist_threshold` (numeric): Defines when a pathway is considered antagonistic. Default: `2.0`.
#'       }
#'     \item **For Variance Filtering (`"variance"`)**:
#'       \itemize{
#'         \item `variance_percentile` (numeric): Percentile threshold for gene variance filtering. Default: `0.2` (removes genes in the lowest 20% of variance).
#'       }
#'     \item **For Correlation Filtering (`"correlation"`)**:
#'       \itemize{
#'         \item `correlation_threshold` (numeric): Proportion of negatively correlated genes required to exclude a pathway. Default: `0.5`.
#'         \item `negative_corr_cutoff` (numeric): Correlation coefficient threshold below which genes are considered antagonistic. Default: `-0.3`.
#'       }
#'   }
#' @param log_transform Logical. Whether to log10 transform the dose values. Default is `FALSE`.
#' @param spline_knots An integer specifying the number of internal knots for the spline function. Default is `2`.
#'                     The number of knots should be lower than the number of unique dose values.
#' @param knot_method A string specifying the method for selecting spline knots. Options:
#'   \itemize{
#'     \item `"quantile"` (default): Selects knots based on percentiles of the dose distribution.
#'     \item `"geometric"`: Places knots at logarithmically spaced intervals.
#'     \item `"manual"`: Uses user-defined knots from the `manual_knots` parameter.
#'   }
#' @param manual_knots A numeric vector specifying custom knot positions (used only if `knot_method = "manual"`).
#'                     If `NULL`, this option is ignored.
#' @param num_cores The number of cores to use for parallel processing (default is 5).
#'
#' @return A list containing the results of the DoseRider analysis for each gene set.
#' @import SummarizedExperiment
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import progress
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar
#' @importFrom doSNOW registerDoSNOW
#'
#' @examples
#' \dontrun{
#' data("SummarizedExperiment")
#' gmt <- list(geneSet1 = list(genes = c("gene1", "gene2")))
#' results <- DoseRiderParallel(se, gmt, "dose", "sample", "covariate", "rnaseq", num_cores = 4, modelType = "GAMM")
#' }
#'
#' @export
DoseRiderParallel <- function(se, gmt, dose_col = "dose", sample_col = "sample",
                              covariates = c(), omic = "rnaseq", minGSsize = 5,
                              maxGSsize = 300, method = "fdr", num_cores = 5,
                              modelType = "LMM",
                              FilterPathway = T,
                              filters = c("correlation"),
                              filter_params = list(),
                              log_transform = F,
                              models = c("linear", "non_linear_fixed","non_linear_mixed"),
                              spline_knots = 2, knot_method = "quantile",
                              manual_knots = NULL) {
  # Register the parallel backend
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)

  # Validate input data and metadata columns
  validate_input_doserider(se, dose_col, sample_col, covariates)

  # Initialize progress bar and options for parallel processing
  total_gene_sets <- length(gmt)
  pb <- txtProgressBar(min = 0, max = total_gene_sets, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))

  # Loop over gene sets in parallel
  results <- foreach(i = seq_along(gmt), .packages = c("SummarizedExperiment", "lme4", "doseRider", "dplyr"),
                     .combine = 'c', .options.snow = opts) %dopar% {
                       geneset_results <- suppressWarnings(process_gene_set(
                         se = se, dose_col = dose_col, sample_col = sample_col, omic = omic, gmt = gmt,
                         i = i, minGSsize = minGSsize, maxGSsize = maxGSsize, covariates = covariates,
                         modelType = modelType,
                         FilterPathway = FilterPathway, filters = filters,
                         filter_params = filter_params,
                         log_transform = log_transform, models = models,
                         spline_knots = spline_knots, knot_method = knot_method,
                         manual_knots = manual_knots
                       ))
                       if (!is.null(geneset_results)) {
                         setNames(list(geneset_results), gmt[[i]]$pathway)
                       } else {
                         NULL
                       }
                     }

  # Close parallel backend and progress bar
  stopCluster(cl)
  close(pb)

  # Adjust p-values for multiple testing
  results <- adjust_pvalues_doserider_result(results = results, method = method)

  class(results) <- "DoseRider"
  return(results)
}


#' Filter DoseRider Object
#'
#' This function filters a DoseRider object based on the specified model type and p-value criteria.
#'
#' @param doseRiderObj A DoseRider object containing results from dose-response analysis.
#' @param model_type The type of model to filter by. Options are "all", "linear", "non_linear" (which includes "non_linear_fixed" and "non_linear_mixed").
#' @param filter_type The type of p-value to filter by. Options are "pvalue" for raw p-value and "fdr" for adjusted p-value.
#' @param threshold The p-value threshold for filtering. Defaults to 0.05.
#'
#' @return A filtered DoseRider object.
#'
#' @examples
#' \dontrun{
#' filtered_results <- filter_DoseRider(doseRiderObj, model_type = "linear", filter_type = "fdr", threshold = 0.05)
#' }
#'
#' @export
filter_DoseRider <- function(doseRiderObj, model_type = "all", filter_type = "pvalue", threshold = 0.05) {
  valid_model_types <- c("all", "linear", "non_linear", "non_linear_fixed", "non_linear_mixed")
  valid_filter_types <- c("pvalue", "fdr")

  # Check for valid model_type and filter_type
  if (!model_type %in% valid_model_types) {
    stop("Invalid model_type. Choose from 'all', 'linear', 'non_linear', 'non_linear_fixed', 'non_linear_mixed'.")
  }

  if (!filter_type %in% valid_filter_types) {
    stop("Invalid filter_type. Choose from 'pvalue' or 'fdr'.")
  }

  # Function to determine if a result meets the filter criteria
  meets_criteria <- function(result) {
    # Select the p-value column based on filter_type
    pvalue_col <- switch(filter_type,
                         "pvalue" = "best_model_pvalue",
                         "fdr" = "best_model_adj_pvalue")

    # Check if the result meets the p-value threshold
    if (!is.null(result[[pvalue_col]]) && !is.na(result[[pvalue_col]]) && result[[pvalue_col]] <= threshold) {

      # Filtering logic based on model_type
      if (model_type == "all") {
        return(TRUE)
      } else if (model_type == "linear" && result$best_model == "linear") {
        return(TRUE)
      } else if (model_type == "non_linear" &&
                 (result$best_model == "non_linear_fixed" || result$best_model == "non_linear_mixed")) {
        return(TRUE)
      } else if (model_type == result$best_model) {  # for specific models like "non_linear_fixed" or "non_linear_mixed"
        return(TRUE)
      }
    }

    return(FALSE)
  }

  # Filter the DoseRider object
  filtered_results <- lapply(doseRiderObj, function(result) {
    if (meets_criteria(result)) {
      return(result)
    } else {
      return(NULL)
    }
  })

  # Remove NULL entries
  filtered_results <- Filter(Negate(is.null), filtered_results)

  # Return as a DoseRider object
  class(filtered_results) <- "DoseRider"
  return(filtered_results)
}




#' Convert a DoseRider Object to a Data Frame
#'
#' This function converts a DoseRider object into a data frame for easier analysis and visualization.
#' It extracts various attributes associated with dose-response analysis results, handling nested list structures.
#'
#' @param object A DoseRider object containing results from dose-response analysis.
#'
#' @return A data frame with attributes from the DoseRider object, where each row corresponds to one gene set.
#'
#' @examples
#' \dontrun{
#'   # Assuming `dose_rider_result` is a DoseRider object
#'   result_df <- as.data.frame.DoseRider(dose_rider_result)
#' }
#'
#' @export
as.data.frame.DoseRider <- function(object) {
  # Define the attributes to extract
  # Define the attributes to extract
  attrs_to_extract <- c("Geneset","Geneset_Size","Genes",
                        "null_AIC", "null_AICc","null_BIC","null_df","null_ICC",
                        "linear_AIC","linear_AICc", "linear_BIC", "linear_df","linear_ICC",
                        "non_linear_fixed_AIC","non_linear_fixed_AICc","non_linear_fixed_BIC","non_linear_fixed_df", "non_linear_fixed_ICC",
                        "non_linear_mixed_AIC", "non_linear_mixed_AICc","non_linear_mixed_BIC","non_linear_mixed_df","non_linear_mixed_ICC",
                        "p_value_linear","p_value_non_linear_fixed","p_value_non_linear_mixed",
                        "adjusted_linear_p_value","adjusted_non_linear_fixed_p_value","adjusted_non_linear_mixed_p_value",
                        "best_model", "best_model_pvalue", "best_model_adj_pvalue","best_model_aicc", "best_model_bic",
                        "OptimalClusters" )


  # Convert the DoseRider object to a list of data frames, handling nested lists
  df_list <- lapply(object, function(x) {
    converted_df <- convert_nested_list_to_df(x[attrs_to_extract])
    return(converted_df)
  })

  # Combine individual data frames into one
  results_df <- do.call(rbind, df_list)
  results_df <- results_df[!is.na(results_df$Geneset),]
  return(results_df)
}

# Utility Function: Convert Nested List to DataFrame
convert_nested_list_to_df <- function(lst) {
  # Find the maximum length among the list elements
  max_len <- max(sapply(lst, length))

  # Ensure each list element has the same length by padding shorter elements with NA
  lst_padded <- lapply(lst, function(el) {
    length_diff <- max_len - length(el)
    if (length_diff > 0) {
      el <- c(el, rep(NA, length_diff))
    }
    return(el)
  })

  # Convert the padded list to a data frame
  df <- as.data.frame(lst_padded, stringsAsFactors = FALSE)
  return(df)
}



#' Print method for DoseRider
#'
#' This function prints a summary of a DoseRider object.
#' It shows the class of the object, the number of gene sets it contains, and some details about the first few gene sets.
#'
#' @param x A DoseRider object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#' \dontrun{
#'   # Assuming `dose_rider_result` is a DoseRider object
#'   print(dose_rider_result)
#' }
#'
#' @export
print.DoseRider <- function(x, ...) {
  cat("Object of class DoseRider\n")
  cat(paste("Number of gene sets:", length(x), "\n"))

  if (length(x) > 0) {
    cat("Summary of the first few gene sets:\n")
    max_display <- min(3, length(x))
    for (i in seq_len(max_display)) {
      geneset_name <- names(x)[i]
      geneset_info <- x[[geneset_name]]
      cat(paste("\nGene Set", i, ":", geneset_name, "\n"))
      cat(paste("  Number of genes:", geneset_info$Genes, "\n"))
      cat(paste("  Best model type:", geneset_info$best_model, "\n"))
      cat(paste("  Best model AICc:", geneset_info$best_model_aicc, "\n"))
      cat(paste("  P-value of best model:", geneset_info$best_model_pvalue, "\n"))
      if (!is.null(geneset_info$OptimalClusters)) {
        cat(paste("  Optimal number of clusters:", geneset_info$OptimalClusters, "\n"))
      }
    }

    if (length(x) > max_display) {
      cat(paste("\n...and", length(x) - max_display, "more gene sets.\n"))
    }
  }
}

