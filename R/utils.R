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
  return(results)
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



