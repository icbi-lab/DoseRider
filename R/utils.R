# Validate input data
validate_input_doserider <- function(se,dose_col,sample_col,covariates){
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input 'se' should be a SummarizedExperiment object.")
  }
  metadata <- colData(se)
  required_columns <- c(dose_col, sample_col, covariates)
  missing_columns <- setdiff(required_columns, names(metadata))
  if (length(missing_columns) > 0) {
    stop("Missing columns in metadata: ", paste(missing_columns, collapse = ", "))
  }
}


# Adjust p-values for multiple testing
adjust_pvalues_doserider_result <- function(results, method){
  non_linear_fixed_p_values <- sapply(results, function(x) x$p_value_non_linear_fixed)
  non_linear_mixed_p_values <- sapply(results, function(x) x$p_value_non_linear_mixed)
  linear_p_values <- sapply(results, function(x) x$p_value_linear)

  adjusted_non_linear_fixed_p_values <- p.adjust(non_linear_fixed_p_values, method)
  adjusted_non_linear_mixed_p_values <- p.adjust(non_linear_mixed_p_values, method)
  adjusted_linear_p_values <- p.adjust(linear_p_values, method)

  # Add adjusted p-values to results
  for (i in seq_along(results)) {
    results[[i]]$adjusted_non_linear_fixed_p_value <- adjusted_non_linear_fixed_p_values[i]
    results[[i]]$adjusted_non_linear_mixed_p_value <- adjusted_non_linear_mixed_p_values[i]
    results[[i]]$adjusted_linear_p_value <- adjusted_linear_p_values[i]
  }
  # Add best_model_padj coulmn
  results <- add_best_model_adj_pvalue(results)
  return(results)
}

#' Add Best Model Adjusted P-value to DoseRider Object
#'
#' This function adds a column "best_model_adj_pvalue" to a DoseRider object, containing the adjusted p-value
#' corresponding to the best model for each gene set.
#'
#' @param doseRiderObj A DoseRider object containing results from dose-response analysis.
#'
#' @return A DoseRider object with an added column "best_model_adj_pvalue".
#'
#' @examples
#' \dontrun{
#' updated_results <- add_best_model_adj_pvalue(doseRiderObj)
#' }
#'
#' @export
add_best_model_adj_pvalue <- function(doseRiderObj) {

  for (i in seq_along(doseRiderObj)) {
    result <- doseRiderObj[[i]]
    best_model <- result$best_model

    if (best_model == "linear") {
      result$best_model_adj_pvalue <- result$adjusted_linear_p_value
    } else if (best_model == "non_linear_fixed") {
      result$best_model_adj_pvalue <- result$adjusted_non_linear_fixed_p_value
    } else if (best_model == "non_linear_mixed") {
      result$best_model_adj_pvalue <- result$adjusted_non_linear_mixed_p_value
    } else {
      result$best_model_adj_pvalue <- NA
    }

    doseRiderObj[[i]] <- result
  }

  class(doseRiderObj) <- "DoseRider"
  return(doseRiderObj)
}

#' Function to obtain the top k gene sets sorted by a specified column
#'
#' @param dose_rider_results Data frame containing gene set information
#' @param pvalue_column Column name for p-values (default: "best_model_adj_pvalue")
#' @param order_column Column name for ordering (default: "NegLogPValue")
#' @param top Number of top gene sets to retrieve (default: 10)
#' @param decreasing Logical indicating sort order (default: TRUE)
#' @return Vector of gene set names of the top k gene sets

get_top_genesets <- function(dose_rider_results,
                             pvalue_column = "best_model_adj_pvalue",
                             order_column = "NegLogPValue",
                             top = 10,
                             decreasing = TRUE) {
  # Convert to data frame if needed
  dose_rider_df <- as.data.frame.DoseRider(dose_rider_results)

  # Add -log10(p-value) column
  dose_rider_df[["NegLogPValue"]] <- -log10(dose_rider_df[[pvalue_column]])

  # Sort by the specified column in the desired order
  dose_rider_df <- dose_rider_df[order(dose_rider_df[[order_column]], decreasing = decreasing), ]

  # Filter out rows with NA values in the ordering column
  dose_rider_df <- dose_rider_df[!is.na(dose_rider_df[[order_column]]), ]

  # Retrieve the top pathways
  top_pathways_df <- head(dose_rider_df, top)

  # Extract gene set names
  gene_set_names <- top_pathways_df$Geneset

  return(gene_set_names)
}


