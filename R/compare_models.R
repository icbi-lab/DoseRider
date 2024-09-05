#' Compare Models Individually Against the Null Model with ANOVA
#'
#' This function compares the goodness of fit of individual models (linear,
#' non-linear fixed effects, and non-linear mixed effects) against a null/base model
#' using ANOVA. P-values are calculated for each model comparison individually.
#'
#' @param null_results A fitted model representing the null or base model.
#' @param linear_results A fitted model representing the linear model (optional).
#' @param non_linear_fixed_results A fitted model representing the non-linear fixed effects model (optional).
#' @param non_linear_mixed_results A fitted model representing the non-linear mixed effects model (optional).
#' @param modelType A string indicating the type of model comparison: "LMM" for linear mixed models or "GAMM" for generalized additive mixed models.
#'
#' @return A named list of p-values for the comparisons of each model against the base model (null_results).
#' If a model is not provided, its p-value is returned as NA.
#'
#' @examples
#' # Assuming lmm_null, lmm_linear, lmm_non_linear_fixed, and lmm_non_linear_mixed
#' # are fitted models:
#' # compare_models_pvalues <- compare_all_models(lmm_null, lmm_linear, lmm_non_linear_fixed,
#' # lmm_non_linear_mixed, modelType = "LMM")
#' @importFrom mgcv anova.gam
compare_all_models <- function(null_results, linear_results = NULL,
                               non_linear_fixed_results = NULL,
                               non_linear_mixed_results = NULL,
                               modelType) {

  # Initialize an empty list to store p-values
  p_values <- list()

  if (modelType == "LMM") {
    # Perform individual ANOVA comparisons for LMM models

    # Compare null vs. linear
    if (!is.null(linear_results)) {
      anova_linear <- anova(null_results, linear_results)
      p_values[["p_value_linear"]] <- anova_linear[2, "Pr(>Chisq)"]
    } else {
      p_values[["p_value_linear"]] <- NA
    }

    # Compare null vs. non-linear fixed
    if (!is.null(non_linear_fixed_results)) {
      anova_non_linear_fixed <- anova(null_results, non_linear_fixed_results)
      p_values[["p_value_non_linear_fixed"]] <- anova_non_linear_fixed[2, "Pr(>Chisq)"]
    } else {
      p_values[["p_value_non_linear_fixed"]] <- NA
    }

    # Compare null vs. non-linear mixed
    if (!is.null(non_linear_mixed_results)) {
      anova_non_linear_mixed <- anova(null_results, non_linear_mixed_results)
      p_values[["p_value_non_linear_mixed"]] <- anova_non_linear_mixed[2, "Pr(>Chisq)"]
    } else {
      p_values[["p_value_non_linear_mixed"]] <- NA
    }

  } else if (modelType == "GAMM") {
    # Perform individual ANOVA comparisons for GAMM models

    # Compare null vs. linear
    if (!is.null(linear_results)) {
      anova_linear <- anova.gam(null_results, linear_results, test = "F")
      p_values[["p_value_linear"]] <- anova_linear[2, "Pr(>F)"]
    } else {
      p_values[["p_value_linear"]] <- NA
    }

    # Compare null vs. non-linear fixed
    if (!is.null(non_linear_fixed_results)) {
      anova_non_linear_fixed <- anova.gam(null_results, non_linear_fixed_results, test = "F")
      p_values[["p_value_non_linear_fixed"]] <- anova_non_linear_fixed[2, "Pr(>F)"]
    } else {
      p_values[["p_value_non_linear_fixed"]] <- NA
    }

    # Compare null vs. non-linear mixed
    if (!is.null(non_linear_mixed_results)) {
      anova_non_linear_mixed <- anova.gam(null_results, non_linear_mixed_results, test = "F")
      p_values[["p_value_non_linear_mixed"]] <- anova_non_linear_mixed[2, "Pr(>F)"]
    } else {
      p_values[["p_value_non_linear_mixed"]] <- NA
    }
  } else {
    stop("Unsupported model type. Choose 'LMM' or 'GAMM'.")
  }

  # Return the list of p-values
  return(p_values)
}


#' Select Best Model Based on AICc and P-Values
#'
#' This function selects the best model from a list of models (null, linear,
#' non-linear fixed, and non-linear mixed) by comparing their AICc values and
#' significance of p-values against the null model. The best model is the one with
#' the lowest AICc and a significant p-value.
#'
#' @param model_list A named list of model objects. Names should include 'null', 'linear',
#' 'non_linear_fixed', and 'non_linear_mixed'.
#' @param p_value_list A named list of p-values for model comparisons against the null model.
#' @param alpha Significance level for p-value comparison (default is 0.05).
#' @return The name of the best fitting model ('null', 'linear', 'non_linear_fixed', or 'non_linear_mixed').
#' @importFrom AICcmodavg AICc
#' @examples
#' \dontrun{
#' # Assuming you have a list of models called 'model_list'
#' # and a list of p-values called 'p_value_list':
#' best_model <- select_best_model(model_list, p_value_list)
#' }
select_best_model <- function(model_list, p_value_list, alpha = 0.05) {

  # Initialize an empty vector to store model names, p-values, and AICc values
  valid_models <- list()

  # Check each model (linear, non-linear fixed, non-linear mixed) against the null
  for (model_name in c("linear", "non_linear_fixed", "non_linear_mixed")) {

    # Check if model is provided and its p-value is significant
    if (!is.null(model_list[[model_name]]) && p_value_list[[paste0("p_value_", model_name)]] < alpha) {

      # Add the model to the valid models list with its AICc value
      valid_models[[model_name]] <- AICc(model_list[[model_name]])
    }
  }

  # Start by assuming the null model is the best
  best_model <- "null"
  best_AICc <- AICc(model_list$null)

  # Find the model with the lowest AICc among valid models
  if (length(valid_models) > 0) {
    # Get the model with the lowest AICc
    lowest_AICc_model <- names(which.min(valid_models))
    lowest_AICc <- min(unlist(valid_models))

    # Check if the lowest AICc is better than the null model (by a delta of at least 2)
    if (lowest_AICc < (best_AICc - 2)) {
      best_model <- lowest_AICc_model
    }
  }

  return(best_model)
}



# The `create_legend_labels` function compiles key model statistics into a concise summary for plot legends from dose-response analysis results.
create_legend_labels <- function(dose_rider_results, gene_set_name) {
  # Extract the best model information
  best_model <- dose_rider_results[[gene_set_name]]$best_model
  p_value <- round(dose_rider_results[[gene_set_name]]$best_model_pvalue, 3)
  p_value_label <- ifelse(p_value < 0.001, "<0.001", as.character(p_value))
  aicc <- round(dose_rider_results[[gene_set_name]]$best_model_aicc, 2)

  # Construct the legend label text
  legend_labels <- paste("Best model:", best_model,
                         ", model p-value:", p_value_label,
                         ", model AICc:", aicc)

  return(legend_labels)
}

