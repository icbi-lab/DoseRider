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


theme_dose_rider <- function(legend_position = "none", text_size = 5, margin_space = 1, fix_ratio = TRUE) {
  # Base theme
  custom_theme <- theme_minimal() +
    theme(
      text = element_text(size = text_size),
      plot.caption = element_text(size = text_size - 2),
      plot.title = element_text(size = text_size, hjust = 0.5),
      axis.title = element_text(size = text_size),
      axis.text = element_text(size = text_size),
      legend.title = element_text(size = text_size),
      legend.text = element_text(size = text_size),
      legend.position = legend_position,
      plot.margin = margin(t = margin_space,  # Top margin
                           r = margin_space,  # Right margin
                           b = margin_space,  # Bottom margin
                           l = margin_space)
    )

  # Conditionally set aspect ratio
  if (fix_ratio) {
    custom_theme <- custom_theme + theme(aspect.ratio = 1)
  }

  return(custom_theme)
}


custom_palette <- c(
  "#D32F2F", "#7B1FA2", "#303F9F", "#0288D1", "#00796B",
  "#388E3C", "#689F38", "#AFB42B", "#FBC02D", "#FFA000",
  "#F57C00", "#E64A19", "#5D4037", "#616161", "#455A64",
  "#C2185B", "#512DA8", "#1976D2", "#0097A7", "#00796B",
  "#8D6E63", "#78909C", "#6D4C41", "#546E7A", "#BF360C",
  "#3E2723", "#0D47A1", "#1B5E20", "#33691E", "#827717",
  "#F57F17", "#FF6F00", "#E65100", "#BF360C", "#3E2723",
  "#263238", "#212121", "#DD2C00", "#FF3D00", "#FF6E40",
  "#FF9E80", "#6A1B9A", "#AB47BC", "#BA68C8", "#CE93D8"
)



