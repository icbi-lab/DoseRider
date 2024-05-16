#' Perform DoseRider Analysis Using Linear Mixed Models (LMMs) or Generalized Additive Mixed Models (GAMMs)
#'
#' This function performs DoseRider analysis on gene expression data, applying either Linear Mixed Models (LMMs) or
#' Generalized Additive Mixed Models (GAMMs) to each gene set defined in the gene matrix transposed (GMT) format.
#' It evaluates dose-response relationships in the context of gene sets and calculates various model metrics,
#' significance, Benchmark Dose (BMD), and smoothing predictions.
#'
#' @param se SummarizedExperiment object or a matrix/data frame containing gene expression data.
#' @param gmt List of gene sets, each represented as a list with gene names.
#' @param dose_col Name of the column representing dose information.
#' @param sample_col Name of the column representing sample information.
#' @param covariates Optional, vector specifying the covariate column(s) in `se`.
#' @param omic Type of omics data, defaults to "rnaseq".
#' @param minGSsize Minimum gene set size for analysis, defaults to 5.
#' @param maxGSsize Maximum gene set size for analysis, defaults to 300.
#' @param method Method for multiple testing adjustment, defaults to "fdr".
#' @param modelType Type of model to fit for each gene set, options are "LMM" for Linear Mixed Models
#'                  and "GAMM" for Generalized Additive Mixed Models. Defaults to "LMM".
#'
#' @return A list containing results for each gene set including various metrics, p-values,
#'         and adjusted p-values. The structure of results will depend on the model type used.
#'
#' @examples
#' \dontrun{
#' data("SummarizedExperiment")
#' gmt <- list(geneSet1 = list(genes = c("gene1", "gene2")))
#' results <- DoseRider(se, gmt, "dose", "sample", "covariate", "rnaseq", modelType = "GAMM")
#' }
#'
#' @importFrom stats p.adjust
#' @import utils

process_gene_set <- function(se, dose_col, sample_col, omic, gmt, i, minGSsize = 5, maxGSsize = 300, covariates = c(), modelType = "LMM") {
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
  # Prepare data for the gene set
  long_df <- suppressWarnings(prepare_data(se=se, geneset=geneset, dose_col=dose_col, sample_col=sample_col, omic=omic))

  # Check gene set size criteria
  if (is.null(long_df) || length(unique(long_df$gene)) < minGSsize || length(unique(long_df$gene)) > maxGSsize) {
    return(NULL)
  }

  # Fit LMM with null and cubic spline models
  #k <-  length(unique(as.vector(long_df[[dose_col]]))) - 1
  if (modelType == "LMM"){
    null_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "null", omic)
    linear_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "linear", omic)
    non_linear_fixed_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "non_linear_fixed", omic)
    non_linear_mixed_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "non_linear_mixed", omic)

    null_results <- suppressWarnings(fit_lmm(null_formula, long_df, omic))
    linear_results <- suppressWarnings(fit_lmm(linear_formula, long_df, omic))
    non_linear_fixed_results <- suppressWarnings(fit_lmm(non_linear_fixed_formula, long_df, omic))
    non_linear_mixed_results <- suppressWarnings(fit_lmm(non_linear_mixed_formula, long_df, omic))

  } else if (modelType == "GAMM"){
    null_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "null", omic)
    linear_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "linear", omic)
    non_linear_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "non_linear", omic)
    non_linear_fixed_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "non_linear_fixed", omic)
    non_linear_mixed_formula <- create_gamm_formula("counts", dose_col, "gene", covariates, "non_linear_mixed", omic)

    null_results <- suppressWarnings(fit_gam(null_formula, long_df, omic))
    linear_results <- suppressWarnings(fit_gam(linear_formula, long_df, omic))
    non_linear_fixed_results <- suppressWarnings(fit_gam(non_linear_fixed_formula, long_df, omic))
    non_linear_mixed_results <- suppressWarnings(fit_gam(non_linear_mixed_formula, long_df, omic))
  }

  # Compare models and compute best model
  fitted_models <- list(null = null_results, linear = linear_results,
                        non_linear_fixed = non_linear_fixed_results,
                        non_linear_mixed = non_linear_mixed_results)
  is_fitted <- sapply(fitted_models, is_fitted_model)
  fitted_models <- fitted_models[is_fitted]

  if (length(fitted_models) > 1) {
    p_value_list <- compare_all_models(null_results, linear_results,
                                       non_linear_fixed_results,
                                       non_linear_mixed_results, modelType)
    best_model_AICc <- select_best_model(fitted_models, p_value_list)
    best_model <- fitted_models[[best_model_AICc]]
  } else {
    p_value_list <- list("p_value_linear"= NA, "p_value_non_linear"= NA)
    best_model_AICc <- names(fitted_models)
    best_model <- fitted_models[[1]]
  }

  # Compute metrics for null and cubic models
  null_metrics <- if (modelType=="LMM") compute_metrics_lmm(null_results) else compute_metrics_gamm(null_results)
  linear_metrics <- if (modelType=="LMM") compute_metrics_lmm(linear_results) else compute_metrics_gamm(linear_results)
  non_linear_fixed_metrics <- if (modelType=="LMM") compute_metrics_lmm(non_linear_fixed_results) else compute_metrics_gamm(non_linear_fixed_results)
  non_linear_mixed_metrics <- if (modelType=="LMM") compute_metrics_lmm(non_linear_mixed_results) else compute_metrics_gamm(non_linear_mixed_results)

  # Calculate smoothing values and BMD for the best model
  if (best_model_AICc != "null") {
    smooth_values <- smooth_pathway_trend(best_model, long_df, dose_col, sample_col, omic, TRUE, covariates, dose_points = 50)
    random_effect <- if (modelType == "LMM") extract_random_effects_lmm(best_model, dose_col) else extract_random_effects_gamm(best_model, dose_col)
    long_df$predictions <- predict(best_model, newdata = long_df)
    optimal_clusters_silhouette <- optimal_clusters_silhouette(smooth_values, dose_col, max_clusters = 10)
    cluster <- optimal_clusters_silhouette$Cluster
    n_cluster <- optimal_clusters_silhouette$OptimalClusters
    cluster_specific_results <- list()

    #Compute mean trend, bmd and derivates for each cluster
    for (j in c(1:n_cluster)) {
      cluster_genes <- names(cluster[cluster == j])
      smooth_cluster <- smooth_values[smooth_values$gene %in% cluster_genes,]
      derivative_cluster <- compute_derivatives(smooth_cluster, dose_col)
      bmd_cluster <- compute_bmd_from_main_trend(smooth_cluster, dose_col, z = 1)
      # Store the cluster-specific results
      # Store the cluster-specific results
      cluster_specific_results[[paste("Cluster", j)]] <- list(
        Derivative = derivative_cluster,
        #SmoothPathway = smooth_pathway_cluster,
        BMD = bmd_cluster
      )
    }

  } else {
    smooth_values <- NA
    smooth_pathway <- NA
    derivative <- NA
    #bmd <- NA
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
    ClusterSpecificResults = cluster_specific_results
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
#' @param se SummarizedExperiment object or a matrix/data frame containing gene expression data.
#' @param gmt List of gene sets, each represented as a list with gene names.
#' @param dose_col Name of the column representing dose information.
#' @param sample_col Name of the column representing sample information.
#' @param covariates Optional, vector specifying the covariate column(s) in `se`.
#' @param omic Type of omics data, defaults to "rnaseq".
#' @param minGSsize Minimum gene set size for analysis, defaults to 5.
#' @param maxGSsize Maximum gene set size for analysis, defaults to 300.
#' @param method Method for multiple testing adjustment, defaults to "fdr".
#' @param modelType Type of model to be used for analysis, "LMM" for Linear Mixed Models or "GAMM"
#' for Generalized Additive Mixed Models. Defaults to "LMM".
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
                         maxGSsize = 300, method = "fdr", modelType = "LMM") {

  # Validate input data
  validate_input_doserider(se,dose_col,sample_col,covariates)

  # Initialize results list
  results <- list()
  total_gene_sets <- length(gmt)
  pb <- txtProgressBar(min = 0, max = total_gene_sets, style = 3)

  # Process each gene set
  for (i in seq_along(gmt)) {
    setTxtProgressBar(pb, i)
    geneset_results <- suppressMessages(suppressWarnings(process_gene_set(se, dose_col,
                                        sample_col, omic, gmt, i, minGSsize, maxGSsize, covariates, modelType)))
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
#' @param se The input SummarizedExperiment object or matrix/data frame with metadata.
#' @param gmt The gene set collection as a list with gene sets.
#' @param dose_col The name of the column in the metadata representing the dose information.
#' @param sample_col The name of the column in the metadata representing the sample information.
#' @param covariates The name of the column in the metadata representing the covariate information (optional).
#' @param omic The type of omic data used (default is "rnaseq").
#' @param minGSsize The minimum gene set size for filtering (default is 5).
#' @param maxGSsize The maximum gene set size for filtering (default is 300).
#' @param method The p-value adjustment method for FDR correction (default is "fdr").
#' @param modelType Type of model to be used for analysis, "LMM" for Linear Mixed Models or "GAMM"
#' for Generalized Additive Mixed Models. Defaults to "LMM".
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
#'
#' @export
DoseRiderParallel <- function(se, gmt, dose_col = "dose", sample_col = "sample",
                                 covariates = c(), omic = "rnaseq", minGSsize = 5,
                                 maxGSsize = 300, method = "fdr", num_cores = 5,
                                 modelType = "LMM") {
  # Register the parallel backend
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)

  # Validate input data and metadata columns
  validate_input_doserider(se,dose_col,sample_col,covariates)

  # Initialize progress bar and options for parallel processing
  total_gene_sets <- length(gmt)
  pb <- txtProgressBar(min = 0, max = total_gene_sets, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))

  # Loop over gene sets in parallel
  results <- foreach(i = seq_along(gmt), .packages = c("SummarizedExperiment", "lme4", "doseRider","dplyr"),
                     .combine = 'c', .options.snow = opts) %dopar% {
                       geneset_results <- suppressWarnings(process_gene_set(se, dose_col, sample_col,
                                                            omic, gmt, i, minGSsize, maxGSsize, covariates, modelType))
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
  attrs_to_extract <- geneset_results_fields <- c("Geneset","Geneset_Size","Genes",
    "null_AIC", "null_AICc","null_BIC","null_df",
    "linear_AIC","linear_AICc", "inear_BIC", "linear_df","linear_ICC",
    "non_linear_fixed_AIC","non_linear_fixed_AICc","non_linear_fixed_BIC","non_linear_fixed_df", "non_linear_fixed_ICC",
    "non_linear_mixed_AIC", "non_linear_mixed_AICc","non_linear_mixed_BIC","non_linear_mixed_df","non_linear_mixed_ICC",
    "p_value_inear","p_value_non_linear_fixed","p_value_non_linear_mixed",
    "adjusted_linear_p_value","adjusted_non_linear_fixed_p_value","adjusted_non_linear_mixed_p_value",
    "best_model", "best_model_pvalue", "best_model_aicc", "best_model_bic",
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
#' It shows the class of the object and the number of gene sets it contains.
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
}
