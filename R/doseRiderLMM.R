#' Process a Single Gene Set for DoseRider Analysis
#'
#' This function processes a single gene set from a gene set collection for DoseRider analysis. It fits
#' Linear Mixed Models (LMMs) with cubic splines to the gene expression data, calculates model metrics,
#' significance, Benchmark Dose (BMD), and smoothing values. The results are compiled into a list.
#'
#' @param se SummarizedExperiment object containing the gene expression data.
#' @param dose_col Character string specifying the dose column name in `se`.
#' @param sample_col Character string specifying the sample column name in `se`.
#' @param omic Character string indicating the type of omics data.
#' @param gmt List of gene sets (Gene Matrix Transposed format).
#' @param i Index of the current gene set in the `gmt` list.
#' @param minGSsize Minimum size of the gene set for analysis (default 5).
#' @param maxGSsize Maximum size of the gene set for analysis (default 300).
#' @param covariates Vector specifying the covariate column(s) in `se`.
#'
#' @return A list containing analysis results for the gene set, or NULL if the gene set size criteria are not met.
#'
#' @importFrom lme4 ranef
#' @examples
#' \dontrun{
#' # Assuming `se` is a SummarizedExperiment object, `gmt` is a list of gene sets
#' gene_set_results <- process_gene_set_lmm(se, "dose", "sample", "rnaseq", gmt, 1, 5, 300, "sample")
#' }
#'
#' @export
process_gene_set_lmm <- function(se, dose_col, sample_col, omic, gmt, i, minGSsize = 5, maxGSsize = 300, covariates = c()) {
  geneset <- gmt[[i]]$genes
  # Prepare data for the gene set
  long_df <- suppressWarnings(prepare_data(se=se, geneset=geneset, dose_col=dose_col, sample_col=sample_col, omic=omic))

  # Check gene set size criteria
  if (is.null(long_df) || length(unique(long_df$gene)) < minGSsize || length(unique(long_df$gene)) > maxGSsize) {
    return(NULL)
  }

  # Fit LMM with null and cubic spline models
  null_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "null", omic)
  linear_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "linear", omic)
  cubic_formula <- create_lmm_formula("counts", dose_col, "gene", covariates, "cubic", omic, k = length(unique(as.vector(long_df[dose_col]))))

  null_results <- suppressWarnings(fit_lmm(null_formula, long_df, omic))
  linear_results <- suppressWarnings(fit_lmm(linear_formula, long_df, omic))
  cubic_results <- suppressWarnings(fit_lmm(cubic_formula, long_df, omic))

  # Compare models and compute best model
  if (!is.na(null_results) && !is.na(linear_results) && !is.na(cubic_results)) {
    p_value_list <- compare_all_models(null_results, linear_results, cubic_results)
    best_model_AICc <- select_best_model(list(null = null_results, linear = linear_results, cubic = cubic_results))
    random_effect <- ranef(cubic_results)$gene
  } else {
    p_value_list <- NA
    best_model_AICc <- "null"
    random_effect <- NA
  }

  # Compute metrics for null and cubic models
  null_metrics <- compute_metrics_lmm(null_results)
  linear_metrics <- compute_metrics_lmm(linear_results)
  cubic_metrics <- compute_metrics_lmm(cubic_results)

  # Calculate smoothing values and BMD
  if (!is.na(cubic_results) && best_model_AICc == "cubic") {
    smooth_values <- smooth_pathway_trend(cubic_results, long_df, dose_col, sample_col, omic, TRUE, covariates)
    smooth_pathway <- smooth_pathway_trend(cubic_results, long_df, dose_col, sample_col, omic, FALSE, covariates,dose_points = 50)
    derivate <- compute_derivatives(cubic_results, long_df, dose_col, omic, covariates)
    bmd <- compute_bmd_from_main_trend(cubic_results, long_df, dose_col, omic, covariates = covariates)
  } else {
    smooth_values <- NA
    smooth_pathway <- NA
    derivate <- NA
    bmd <- NA
  }


  # Compile results into a list
  geneset_results <- list(
    Geneset = gmt[[i]]$pathway,
    Geneset_Size = length(geneset),
    Genes = length(unique(long_df$gene)),
    Null_AIC = null_metrics$AIC,
    Null_AICc = null_metrics$AICc,
    Null_BIC = null_metrics$BIC,
    Null_df = null_metrics$edf,
    Linear_AIC = linear_metrics$AIC,       # Adding Linear model AIC
    Linear_AICc = linear_metrics$AICc,     # Adding Linear model AICc
    Linear_BIC = linear_metrics$BIC,       # Adding Linear model BIC
    Linear_df = linear_metrics$edf,        # Adding Linear model degrees of freedom
    Cubic_AIC = cubic_metrics$AIC,         # Adding Cubic model AIC
    Cubic_AICc = cubic_metrics$AICc,       # Adding Cubic model AICc
    Cubic_BIC = cubic_metrics$BIC,         # Adding Cubic model BIC
    Cubic_df = cubic_metrics$edf,          # Adding Cubic model degrees of freedom
    P_Value_Linear = p_value_list[1],
    P_Value_Cubic = p_value_list[2],
    Best_Model_AICc = best_model_AICc,
    Smooth_Predictions = list(smooth_values),
    Smooth_Predictions_Pathway = list(smooth_pathway),
    Raw_Values = list(long_df),
    BMD = bmd,
    TCD = derivate,
    random_effect = random_effect
  )

  return(geneset_results)
}


#' Perform DoseRider Analysis Using Linear Mixed Models
#'
#' This function performs DoseRider analysis on gene expression data, applying Linear Mixed Models (LMMs)
#' to each gene set defined in the gene matrix transposed (GMT) format. It evaluates dose-response relationships
#' in the context of gene sets and calculates various model metrics, significance, and smoothing predictions.
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
#'
#' @return A list containing results for each gene set including various metrics, p-values, and adjusted p-values.
#'
#' @examples
#' \dontrun{
#' data("SummarizedExperiment")
#' gmt <- list(geneSet1 = list(genes = c("gene1", "gene2")))
#' results <- DoseRiderLMM(se, gmt, "dose", "sample", "covariate", "rnaseq")
#' }
#'
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar
#' @export
DoseRiderLMM <- function(se, gmt, dose_col = "dose", sample_col = "sample",
                         covariates = c(), omic = "rnaseq", minGSsize = 5,
                         maxGSsize = 300, method = "fdr") {
  # Validate input data
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input 'se' should be a SummarizedExperiment object.")
  }

  # Validate columns in the metadata
  validate_columns <- function(metadata, columns) {
    missing_cols <- setdiff(columns, names(metadata))
    if (length(missing_cols) > 0) {
      stop("Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
    }
  }

  metadata <- colData(se)
  validate_columns(metadata, c(dose_col, sample_col))

  if (length(covariates) && !covariates %in% names(metadata)) {
    validate_columns(metadata, covariates)
  }

  # Initialize results list
  results <- list()
  total_gene_sets <- length(gmt)
  pb <- txtProgressBar(min = 0, max = total_gene_sets, style = 3)

  # Process each gene set
  for (i in seq_along(gmt)) {
    setTxtProgressBar(pb, i)
    geneset_results <- suppressMessages(suppressWarnings(process_gene_set_lmm(se, dose_col, sample_col, omic, gmt, i, minGSsize, maxGSsize, covariates)))
    if (!is.null(geneset_results)) {
      results[[gmt[[i]]$pathway]] <- geneset_results
    }
  }

  close(pb)

  # Adjust p-values for multiple testing
  cubic_p_values <- sapply(results, function(x) x$P_Value_Cubic)
  linear_p_values <- sapply(results, function(x) x$P_Value_Linear)

  adjusted_cubic_p_values <- p.adjust(cubic_p_values, method)
  adjusted_linear_p_values <- p.adjust(linear_p_values, method)

  # Add adjusted p-values to results
  for (i in seq_along(results)) {
    results[[i]]$Adjusted_Cubic_P_Value <- adjusted_cubic_p_values[i]
    results[[i]]$Adjusted_Linear_P_Value <- adjusted_linear_p_values[i]
  }

  class(results) <- "DoseRider"
  return(results)
}

#' DoseRiderParallelLMM Function
#'
#' Perform DoseRider analysis in parallel using Linear Mixed Models and multiple cores.
#'
#' @param se SummarizedExperiment object or a matrix/data frame with metadata.
#' @param gmt List of gene sets in Gene Matrix Transposed format.
#' @param dose_col Name of the column representing dose information in metadata.
#' @param sample_col Name of the column representing sample information in metadata.
#' @param covariates Vector of covariate column names in metadata (optional).
#' @param omic Type of omic data (default: "rnaseq").
#' @param minGSsize Minimum gene set size for inclusion (default: 5).
#' @param maxGSsize Maximum gene set size for inclusion (default: 300).
#' @param method Method for p-value adjustment (default: "fdr").
#' @param num_cores Number of cores for parallel processing (default: 5).
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
#'   results <- DoseRiderParallelLMM(se, gmt, dose_col = "dose", sample_col = "sample",
#'                                   covariate = c("age"), num_cores = 8)
#' }
#' @export
DoseRiderParallelLMM <- function(se, gmt, dose_col = "dose", sample_col = "sample",
                                 covariates = c(), omic = "rnaseq", minGSsize = 5,
                                 maxGSsize = 300, method = "fdr", num_cores = 5) {
  # Register the parallel backend
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)

  # Validate input data and metadata columns
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input 'se' should be a SummarizedExperiment object.")
  }
  metadata <- colData(se)
  required_columns <- c(dose_col, sample_col, covariates)
  missing_columns <- setdiff(required_columns, names(metadata))
  if (length(missing_columns) > 0) {
    stop("Missing columns in metadata: ", paste(missing_columns, collapse = ", "))
  }

  # Initialize progress bar and options for parallel processing
  total_gene_sets <- length(gmt)
  pb <- txtProgressBar(min = 0, max = total_gene_sets, style = 3)
  opts <- list(progress = function(n) setTxtProgressBar(pb, n))

  # Loop over gene sets in parallel
  results <- foreach(i = seq_along(gmt), .packages = c("SummarizedExperiment", "lme4", "doseRider"),
                     .combine = 'c', .options.snow = opts) %dopar% {
                       geneset_results <- suppressWarnings(process_gene_set_lmm(se, dose_col, sample_col, omic, gmt, i, minGSsize, maxGSsize, covariates))
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
  cubic_p_values <- sapply(results, function(x) x$P_Value_Cubic)
  linear_p_values <- sapply(results, function(x) x$P_Value_Linear)

  adjusted_cubic_p_values <- p.adjust(cubic_p_values, method)
  adjusted_linear_p_values <- p.adjust(linear_p_values, method)

  # Add adjusted p-values to results
  for (i in seq_along(results)) {
    results[[i]]$Adjusted_Cubic_P_Value <- adjusted_cubic_p_values[i]
    results[[i]]$Adjusted_Linear_P_Value <- adjusted_linear_p_values[i]
  }

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
    "Null_AIC", "Null_AICc","Null_BIC","Null_df", "Linear_AIC","Linear_AICc", "Linear_BIC",
    "Linear_df","Cubic_AIC","Cubic_AICc","Cubic_BIC","Cubic_df","P_Value_Linear","P_Value_Cubic",
    "Best_Model_AICc", "Adjusted_Cubic_P_Value","Adjusted_Linear_P_Value" )


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
