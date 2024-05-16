#' Compare Models with ANOVA
#'
#' This function compares the goodness of fit of three models (base, linear, and cubic)
#' using an ANOVA test and returns the p-values associated with the model comparisons.
#'
#' @param base_results A fitted model representing the base model.
#' @param linear_results A fitted model representing the linear model.
#' @param cubic_results A fitted model representing the cubic model.
#'
#' @return A named vector of p-values for the comparisons: base vs. linear and linear vs. cubic.
#'
#' @examples
#' # Assuming lm_base, lm_linear, and lm_cubic are fitted linear models:
#' # compare_models_pvalues <- compare_three_models(lm_base, lm_linear, lm_cubic)
#' @importFrom mgcv anova.gam
compare_all_models <- function(null_results, linear_results, non_linear_fixed_results,
                               non_linear_mixed_results, modelType) {

  if (modelType == "LMM"){
    # Perform an ANOVA comparison between the three models
    anova_results <- anova(null_results, linear_results,
                           non_linear_fixed_results,
                           non_linear_mixed_results)
  } else if (modelType == "GAMM"){
    anova_results <- anova.gam(null_results, linear_results,
                                non_linear_fixed_results,
                                non_linear_mixed_results, test = "F")

    rownames(anova_results) <- c("null_results","linear_results","non_linear_fixed_results",
                                 "non_linear_mixed_results")
  }
  # Extract p-values for the comparisons
  p_value_linear <- anova_results["linear_results",]$`Pr(>Chisq)`
  p_value_non_linear_fixed <- anova_results["non_linear_fixed_results",]$`Pr(>Chisq)`
  p_value_non_linear_mixed <- anova_results["non_linear_mixed_results",]$`Pr(>Chisq)`

  # Return the p-values
  return(list("p_value_linear"=p_value_linear, "p_value_non_linear_fixed"=p_value_non_linear_fixed,
             "p_value_non_linear_mixed"=p_value_non_linear_mixed))
}

#' Select Best Model Based on AICc and P-Values
#'
#' Compares null, linear, and non-linear models sequentially using AICc and p-values.
#' Starts with null vs. linear; if linear is significantly better, compares linear to non-linear.
#' The best model has the lowest AICc and significant p-value improvement over predecessors.
#'
#' @param model_list A list of lme4 model objects.
#' @param p_value_list A named list of p-values for model comparisons.
#' @param alpha Significance level for p-value comparison (default is 0.05).
#' @return The name of the best fitting model ('null', 'linear', or 'non_linear').
#' @importFrom AICcmodavg AICc
#' @examples
#' \dontrun{
#' # Assuming you have a list of lme4 models called 'model_list'
#' # and a list of p-values called 'p_value_list':
#' best_model <- select_best_model(model_list, p_value_list)
#' }
select_best_model <- function(model_list, p_value_list, alpha = 0.05) {
  # Initialize best model as null
  best_model <- "null"
  # Null vs linear model comparison
  if (p_value_list$p_value_linear < alpha && AICc(model_list$linear) < (AICc(model_list$null) - 2)) {
    best_model <- "linear"

    # Linear vs non-linear model comparison
    if (p_value_list$p_value_non_linear_fixed < alpha && AICc(model_list$non_linear_fixed) < (AICc(model_list$null) - 2)) {
      best_model <- "non_linear_fixed"

      if (p_value_list$p_value_non_linear_mixed < alpha && AICc(model_list$non_linear_mixed) < (AICc(model_list$null) - 2)) {
        best_model <- "non_linear_mixed"

      }
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

