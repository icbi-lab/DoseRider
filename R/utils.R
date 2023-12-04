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
  cubic_p_values <- sapply(results, function(x) x$P_Value_Cubic)
  linear_p_values <- sapply(results, function(x) x$P_Value_Linear)

  adjusted_cubic_p_values <- p.adjust(cubic_p_values, method)
  adjusted_linear_p_values <- p.adjust(linear_p_values, method)

  # Add adjusted p-values to results
  for (i in seq_along(results)) {
    results[[i]]$Adjusted_Cubic_P_Value <- adjusted_cubic_p_values[i]
    results[[i]]$Adjusted_Linear_P_Value <- adjusted_linear_p_values[i]
  }
  return(results)
}

custom_palette <- c(
  "#D32F2F", "#7B1FA2", "#303F9F", "#0288D1", "#00796B",
  "#388E3C", "#689F38", "#AFB42B", "#FBC02D", "#FFA000",
  "#F57C00", "#E64A19", "#5D4037", "#616161", "#455A64",
  "#C2185B", "#7B1FA2", "#512DA8", "#303F9F", "#1976D2",
  "#0288D1", "#0097A7", "#00796B", "#388E3C", "#689F38",
  "#AFB42B", "#FBC02D", "#FFA000", "#F57C00", "#E64A19",
  "#5D4037", "#616161", "#455A64", "#D32F2F", "#C2185B",
  "#7B1FA2", "#512DA8", "#1976D2", "#0288D1", "#0097A7"
)
