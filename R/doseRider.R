#' Process a single gene set for DoseRider analysis
#'
#' This function processes a single gene set from the gene set collection for DoseRider analysis.
#' It fits GAM models to the gene expression data and calculates model metrics and significance.
#' The results are then stored in a list.
#'
#' @param long_df The long-format data frame containing the gene expression data for the current gene set.
#' @param dose_col The name of the column representing the dose information in the long_df.
#' @param sample_col The name of the column representing the sample information in the long_df.
#' @param gmt The gene set collection as a list with gene sets.
#' @param i The index of the current gene set in the gmt list.
#' @param minGSsize The minimum gene set size for filtering (default is 5).
#' @param maxGSsize The maximum gene set size for filtering (default is 300).
#'
#' @return A list containing the results of the DoseRider analysis for the current gene set, or NULL if the gene set does not meet the size criteria.
#'
#' @examples
#' # Example usage of process_gene_set
#' gene_set_results <- process_gene_set(long_df, dose_col = "dose", sample_col = "sample", gmt, i = 1, minGSsize = 5, maxGSsize = 300)
#'
process_gene_set <- function(se, dose_col, sample_col,omic, gmt, i, minGSsize, maxGSsize) {
  geneset <- gmt[[i]]$genes
  # Filter gene set by length
  long_df <- prepare_data(se, geneset, dose_col, sample_col, omic)
  if (!is.null(long_df) && (length(unique(long_df$gene)) >= minGSsize && length(unique(long_df$gene)) <= maxGSsize)) {
    ## Check the minum overlapping genes
    # Compute GeneRatio (if needed)
    # input_vector <- unique(long_df$gene)
    # ratios <- compute_ratios(geneset, all_genes, input_vector)

    # Step 5: Compute the models and get the significance
    # Fit GAM models and compute AIC, BIC, and df
    base_formula <- create_gamm_formula(response = "counts",
                                        fixed_effects = dose_col,
                                        random_effects = "gene",
                                        model_type = "base")
    linear_formula <- create_gamm_formula(response = "counts",
                                          fixed_effects = dose_col,
                                          random_effects = "gene",
                                          model_type = "linear")
    cubic_formula <- create_gamm_formula(response = "counts",
                                         fixed_effects = dose_col,
                                         random_effects = "gene",
                                         model_type = "cubic")

    base_results <- fit_gam(base_formula, long_df)
    linear_results <- fit_gam(linear_formula, long_df)
    cubic_results <- fit_gam(cubic_formula, long_df)

    # Check if base_results and linear_results are not NaN before computing p-value
    p_value_linear <- ifelse(is.list(base_results) & is.list(linear_results),
                             as.numeric(compareGAMM(base_results, linear_results)),
                             NA)

    # Check if base_results and cubic_results are not NaN before computing p-value
    p_value_cubic <- ifelse(is.list(linear_results) & is.list(cubic_results),
                            as.numeric(compareGAMM(base_results, cubic_results)),
                            NA)

    # Get metrics for linear and cubic models
    base_metrics <- compute_metrics(base_results)
    linear_metrics <- compute_metrics(linear_results)
    cubic_metrics <- compute_metrics(cubic_results)

    # Step 6: Calculate smoothing values for all X
    smooth_values <- smooth_predictions(model = cubic_results, long_df = long_df,
                                        dose_col = dose_col, sample_col = sample_col,
                                        covariate = covariate)

    # Step 7: Save the results
    # Create a list to store the results for the current gene set
    geneset_results <- list()

    # Add the results to the gene set list
    geneset_results$Geneset <- gmt[[i]]$pathway
    geneset_results$Geneset_Size <- length(geneset)
    #geneset_results$GeneRatio <- ratios$geneRatio
    #geneset_results$bgRatio <- ratios$bgRatio
    geneset_results$Genes <- length(unique(long_df$gene))
    geneset_results$Base_AIC <- base_metrics$AIC
    geneset_results$Base_BIC <- base_metrics$BIC
    geneset_results$Base_edf <- base_metrics$edf
    geneset_results$Linear_AIC <- linear_metrics$AIC
    geneset_results$Linear_BIC <- linear_metrics$BIC
    geneset_results$Linear_edf <- linear_metrics$edf
    geneset_results$Cubic_AIC <- cubic_metrics$AIC
    geneset_results$Cubic_BIC <- cubic_metrics$BIC
    geneset_results$Cubic_edf <- cubic_metrics$edf
    geneset_results$P_Value_Linear <- p_value_linear
    geneset_results$P_Value_Cubic <- p_value_cubic
    geneset_results$Smooth_Predictions <- list(smooth_values)

    # Return the gene set results
    return(geneset_results)
  }

  # Return NULL if the gene set does not meet the size criteria
  return(NULL)
}


#' DoseRider Function
#'
#' This function performs a series of analysis on gene expression data
#' provided as a SummarizedExperiment object or a matrix/data frame with metadata.
#' It uses Generalized Additive Mixed Models (GAMMs) to analyze the relationship
#' between dose-response and gene sets.
#'
#' @param se The gene expression data as either a SummarizedExperiment object or a matrix/data frame.
#' @param gmt A list of gene sets, each gene set is a list with a 'genes' element containing a character vector of gene symbols.
#' @param dose_col The name of the dose column in the metadata.
#' @param sample_col The name of the sample column in the metadata.
#' @param metadata (Optional) A dataframe that holds metadata information. This is required when 'se' is a matrix or data frame.
#' @param covariate (Optional) The name of a covariate column in the metadata.
#' @param omic The omic type, defaults to "rnaseq".
#' @param minGSsize Minimum gene set size for consideration, defaults to 5.
#' @param maxGSsize Maximum gene set size for consideration, defaults to 300.
#' @param method The method for multiple testing adjustment, defaults to "fdr".
#'
#' @return A list with results from the dose-response analysis for each gene set.
#' Each list contains the following items:
#' - `Geneset`: The name of the gene set.
#' - `Geneset_Size`: The number of genes in the gene set.
#' - `GeneRatio`: The proportion of genes from the gene set present in the data.
#' - `bgRatio`: The proportion of genes in the background set.
#' - `Genes`: The number of unique genes in the long_df.
#' - `Base_AIC`: Akaike Information Criterion for the base model.
#' - `Base_BIC`: Bayesian Information Criterion for the base model.
#' - `Base_edf`: Effective degrees of freedom for the base model.
#' - `Linear_AIC`: Akaike Information Criterion for the linear model.
#' - `Linear_BIC`: Bayesian Information Criterion for the linear model.
#' - `Linear_edf`: Effective degrees of freedom for the linear model.
#' - `Cubic_AIC`: Akaike Information Criterion for the cubic model.
#' - `Cubic_BIC`: Bayesian Information Criterion for the cubic model.
#' - `Cubic_edf`: Effective degrees of freedom for the cubic model.
#' - `P_Value_Linear`: P-value for the linear model.
#' - `P_Value_Cubic`: P-value for the cubic model.
#' - `Smooth_Predictions`: A list containing smoothed predictions from the cubic model.
#' - `Adjusted_P_Value`: P-value for the cubic model adjusted for multiple testing.
#'
#' @export
#' @importFrom stats p.adjust
#' @importFrom utils txtProgressBar

# DoseRider function
DoseRider <- function(se, gmt, dose_col = "dose", sample_col = "sample",
                      covariate = "", omic = "rnaseq", minGSsize=5,
                      maxGSsize=300, method = "fdr") {
  ### Step 1: Data Validation and Metadata Check ###

  # Check if se is a SummarizedExperiment or a matrix/data frame with metadata
  if (inherits(se, "SummarizedExperiment")) {
    num_samples <- ncol(se)
    num_variables <- nrow(se)
    metadata <- as.data.frame(colData(se))
    print(paste("Working with a SummarizedExperiment with", num_samples, "samples and", num_variables, "variables"))
    ## If not check is a matrix
  } else if (is.data.frame(se) || is.matrix(se)) {
    if (is.null(metadata)) {
      stop("Metadata is required when using a matrix or data frame as input. Please provide the metadata.")
    }

    if (!identical(colnames(se), rownames(metadata))) {
      stop("The row names in the matrix/data frame and metadata do not match.")
    }
  } else {
    stop("se should be either a SummarizedExperiment object or a matrix/data frame with metadata.")
  }


  # Check if sample_col, dose_col, and covariate are columns in metadata
  if (!sample_col %in% colnames(metadata)) {
    stop(paste("Column", sample_col, "not found in the metadata. Please provide the correct sample_col."))
  }

  if (!dose_col %in% colnames(metadata)) {
    stop(paste("Column", dose_col, "not found in the metadata. Please provide the correct dose_col."))
  }

  if (covariate != "" && !covariate %in% colnames(metadata)) {
    stop(paste("Column", covariate, "not found in the metadata. Please provide the correct covariate."))
  }

  ### Step 2: Main Processing ###
  # Create an empty data frame to store the results
  results <- list()

  #### Get some values
  # Get the total number of gene sets
  total_gene_sets <- length(gmt)

  # Vector of uniques genes in the db
  all_genes <- unique(unlist(lapply(gmt, `[[`, "genes")))
  input_vector <- unique(rownames(se))

  # Initialize the progress bar
  #pb <- progress_bar$new(total = total_gene_sets)
  pb <- txtProgressBar(min = 0, max = total_gene_sets, style = 3)

  ##Process genesets
  # Loop over gene sets
  for (i in seq_along(gmt)) {
    geneset_results <- process_gene_set(se, dose_col, sample_col,omic, gmt, i, minGSsize, maxGSsize)
    if (!is.null(geneset_results)) {
      # Add the gene set results to the list of results
      results[[gmt[[i]]$pathway]] <- geneset_results
    }
    # Update the progress bar
    setTxtProgressBar(pb, i)
  }

  # Close the progress bar
  close(pb)

  # Add FDR adjustment to the results
  # Get the adjusted p-values
  adjusted_cubic_p_values <- p.adjust(unlist(lapply(results, `[[`, "P_Value_Cubic")), method = method)
  adjusted_linear_p_values <- p.adjust(unlist(lapply(results, `[[`, "P_Value_Linear")), method = method)

  # Save the adjusted p-value inside each gene set
  for (i in seq_along(results)) {
    results[[i]]$Adjusted_Cubic_P_Value <- adjusted_p_values[i]
    results[[i]]$Adjusted_Linear_P_Value <- adjusted_linear_p_values[i]

    }

  # Add the adjusted p-values to the results object
  class(results) <- "DoseRider"

  return(results)
}


#' DoseRiderParallel function
#'
#' Perform DoseRider analysis in parallel using multiple cores.
#'
#' @param se The input SummarizedExperiment object or matrix/data frame with metadata.
#' @param gmt The gene set collection as a list with gene sets.
#' @param dose_col The name of the column in the metadata representing the dose information.
#' @param sample_col The name of the column in the metadata representing the sample information.
#' @param covariate The name of the column in the metadata representing the covariate information (optional).
#' @param omic The type of omic data used (default is "rnaseq").
#' @param minGSsize The minimum gene set size for filtering (default is 5).
#' @param maxGSsize The maximum gene set size for filtering (default is 300).
#' @param method The p-value adjustment method for FDR correction (default is "fdr").
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
#'
#' @export
#'
#' @examples
#' # Example usage of DoseRiderParallel
#' results <- DoseRiderParallel(se, gmt, dose_col = "dose", sample_col = "sample", covariate = "age", num_cores = 8)
#'
# DoseRiderParallel function
DoseRiderParallel <- function(se, gmt, dose_col = "dose", sample_col = "sample",
                              covariate = "", omic = "rnaseq", minGSsize = 5,
                              maxGSsize = 300, method = "fdr", num_cores = 5) {
  # Register the parallel backend
  cl <- makeCluster(num_cores)
  registerDoSNOW(cl)

  ### Step 1: Data Validation and Metadata Check ###
  # Check if se is a SummarizedExperiment or a matrix/data frame with metadata
  if (inherits(se, "SummarizedExperiment")) {
    num_samples <- ncol(se)
    num_variables <- nrow(se)
    metadata <- as.data.frame(colData(se))
    print(paste("Working with a SummarizedExperiment with", num_samples, "samples and", num_variables, "variables"))
    ## If not check is a matrix
  } else if (is.data.frame(se) || is.matrix(se)) {
    if (is.null(metadata)) {
      stop("Metadata is required when using a matrix or data frame as input. Please provide the metadata.")
    }

    if (!identical(colnames(se), rownames(metadata))) {
      stop("The row names in the matrix/data frame and metadata do not match.")
    }
  } else {
    stop("se should be either a SummarizedExperiment object or a matrix/data frame with metadata.")
  }


  # Check if sample_col, dose_col, and covariate are columns in metadata
  if (!sample_col %in% colnames(metadata)) {
    stop(paste("Column", sample_col, "not found in the metadata. Please provide the correct sample_col."))
  }

  if (!dose_col %in% colnames(metadata)) {
    stop(paste("Column", dose_col, "not found in the metadata. Please provide the correct dose_col."))
  }

  if (covariate != "" && !covariate %in% colnames(metadata)) {
    stop(paste("Column", covariate, "not found in the metadata. Please provide the correct covariate."))
  }

  ### Step 2: Main Processing ###
  # Create an empty data frame to store the results
  results <- list()

  #### Get some values
  # Get the total number of gene sets
  total_gene_sets <- length(gmt) + 1

  # Vector of uniques genes in the db
  all_genes <- unique(unlist(lapply(gmt, `[[`, "genes")))
  input_vector <- unique(rownames(se))

  # Initialize the progress bar
  pb <- txtProgressBar(max=total_gene_sets, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)

  # Loop over gene sets in parallel
  results <- foreach(i = seq_along(gmt), .packages = c("mgcv", "tidyverse", "SummarizedExperiment", "reshape", "doseRider"), .combine = "c", .options.snow = opts) %dopar% {
    geneset_results <- process_gene_set(se, dose_col, sample_col, omic, gmt, i, minGSsize, maxGSsize)
    if (!is.null(geneset_results)) {
      setNames(list(geneset_results), gmt[[i]]$pathway)
    } else {
      NULL
    }
  }

  # Close the parallel backend
  stopCluster(cl)

  # Close the progress bar
  close(pb)

  #Add FDR adjustment to the results
  # Get the adjusted p-values
  adjusted_cubic_p_values <- p.adjust(unlist(lapply(results, `[[`, "P_Value_Cubic")), method = method)
  adjusted_linear_p_values <- p.adjust(unlist(lapply(results, `[[`, "P_Value_Linear")), method = method)

  # Save the adjusted p-value inside each gene set
  for (i in seq_along(results)) {
    results[[i]]$Adjusted_Cubic_P_Value <- adjusted_cubic_p_values[i]
    results[[i]]$Adjusted_Linear_P_Value <- adjusted_linear_p_values[i]

  }

  # Add the adjusted p-values to the results object
  class(results) <- "DoseRider"

  return(results)
}



#' Convert a DoseRider object to a data frame
#'
#' This function extracts specific attributes from a DoseRider object and
#' structures them in a data frame for easier manipulation and visualization.
#' The attributes included are: Geneset, Geneset_Size, GeneRatio, bgRatio,
#' Genes, Base_AIC, Base_BIC, Base_edf, Linear_AIC, Linear_BIC, Linear_edf,
#' Cubic_AIC, Cubic_BIC, Cubic_edf, P_Value_Linear, P_Value_Cubic, and Adjusted_P_Value.
#'
#' @param object A DoseRider object.
#'
#' @return A data frame with specific attributes from the DoseRider object.
#' Each row of the data frame corresponds to an element of the DoseRider object.
#'
#' @examples
#' \dontrun{
#'   # Assuming `dose_rider_result` is a DoseRider object
#'   result_df <- as.data.frame.DoseRider(dose_rider_result)
#' }
#'
#' @export
as.data.frame.DoseRider <- function(object) {
  # Define the variables to keep
  keep <- c("Geneset", "Geneset_Size", "GeneRatio", "bgRatio", "Genes", "Base_AIC",
            "Base_BIC", "Base_edf", "Linear_AIC", "Linear_BIC", "Linear_edf",
            "Cubic_AIC", "Cubic_BIC", "Cubic_edf", "P_Value_Linear", "P_Value_Cubic",
            "Adjusted_P_Value" )

  # Apply a function to each element of the list
  df_list <- lapply(object, function(x) x[keep])

  # Convert the list of data frames to a single data frame
  results_df <- as.data.frame(do.call(rbind, df_list))

  return(results_df)
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
