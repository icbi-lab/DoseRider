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
#' @examples
#' # Load example SummarizedExperiment and gene set list
#' data(se)
#' data(gmt)
#'
#' # Run DoseRider
#' results <- DoseRider(se = se, gmt = gmt, dose_col = "dose", sample_col = "sample")
#'
#' # Print the results for a specific gene set
#' print(results[["gene_set_1"]])
#' @importFrom  progress progress_bar
#' @export


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


  ### Step 3: Obtaining the desing formulas for the analysis ###
  #### Obtaining the desing formulas for the analysis
  # Define the formulas for the models
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

  ### Step 4: Loop Over Gene Sets ###
  # filter gmt due size
  # Loop over the gene sets
  for (i in seq_along(gmt)) {
    geneset <- gmt[[i]]$genes
    #print(gmt[[i]]$pathway)

    # Filter gene set by length
    long_df <- prepare_data(se, geneset, dose_col, sample_col, omic)

    ## Check the minum overlapping genes
    if (!is.null(long_df) && (length(unique(long_df$gene)) >= minGSsize && length(unique(long_df$gene)) <= maxGSsize)) {
      #Compute GeneRatio
      #input_vector <- unique(long_df$gene)
      #ratios <- compute_ratios(geneset, all_genes, input_vector)

      ## Step 5: compute the models and get the significance###
      # Fit GAM models and compute AIC, BIC, and df
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

      ### Step 6: calculate smoothing values for all X ###
      smooth_values <- smooth_predictions(model = cubic_results,long_df =long_df,
                                          dose_col = dose_col, sample_col =sample_col,
                                          covariate = covariate)


      ### Step 7: Save the results###
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

      # Add the gene set results to the list of results
      results[[gmt[[i]]$pathway]] <- geneset_results
    }
    # Update the progress bar
    #pb$tick()
    setTxtProgressBar(pb, i)

  }

  # Close the progress bar
  close(pb)

  # Add FDR adjustment to the results
  # Get the adjusted p-values
  adjusted_cubic_p_values <- unlist(lapply(results, `[[`, "P_Value_Cubic"))
  adjusted_linear_p_values <- unlist(lapply(results, `[[`, "P_Value_Linear"))
  
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

#' Summary method for DoseRider
#'
#' This function provides a summary of a DoseRider object, including the total number
#' of gene sets, the number of significant gene sets for linear and cubic splines, and other relevant information.
#'
#' @param object A DoseRider object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @examples
#' \dontrun{
#'   # Assuming `dose_rider_result` is a DoseRider object
#'   summary(dose_rider_result)
#' }

#' @export
summary.DoseRider <- function(object, ...) {
  num_gene_sets <- length(object)
  num_significant_linear <- sum(object$P_Value_Linear < 0.05)
  num_significant_cubic <- sum(object$P_Value_Cubic < 0.05)
  #num_significant_both <- sum(object$Adjusted_P_Value < 0.05)
  
  cat("Summary of DoseRider object:\n")
  cat("=============================\n")
  cat("Total number of gene sets:", num_gene_sets, "\n")
  cat("Number of significant gene sets (Linear):", num_significant_linear, "\n")
  cat("Number of significant gene sets (Cubic):", num_significant_cubic, "\n")
  #cat("Number of significant gene sets (Both):", num_significant_both, "\n")
}

