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
compare_all_models <- function(null_results, linear_results, cubic_results, modelType) {

  if (modelType == "LMM"){
    # Perform an ANOVA comparison between the three models
    anova_results <- anova(null_results, linear_results, cubic_results)
  } else if (modelType == "GAMM"){
    anova_results <- anova.gam(null_results, linear_results,cubic_results, test = "F")
  }
  # Extract p-values for the comparisons
  p_linear_vs_null <- anova_results$Pr[2]
  p_cubic_vs_linear <- anova_results$Pr[3]

  # Return the p-values
  return(c(p_linear_vs_null, p_cubic_vs_linear))
}

#' Select the Best Pathway Model Based on AICc
#'
#' This function evaluates the goodness-of-fit for models that incorporate random effects
#' for the genes within pathways. The function uses the AICc criterion to determine the best
#' model among the provided options. The null model contains only a random intercept for the genes.
#'
#' The best model is chosen as the one giving the lowest AICc value.
#' Models with an AICc value not lower than the AICc of the null model minus 2 are eliminated.
#'
#' @param models A named list of models, where:
#'   - 'null' corresponds to the model with only a random intercept.
#'   - 'base', 'linear', and 'cubic' are models with varying fixed effects, but all incorporate
#'     random effects for the genes in the pathways.
#'
#' @return A character string indicating the name of the best model for the pathway.
#' @importFrom AICcmodavg AICc
#' @examples
#' \dontrun{
#' library(lme4)
#' data <- data.frame(y = rnorm(100), x = rnorm(100), gene = sample(letters[1:5], 100, replace = TRUE))
#' null_model <- lmer(y ~ 1 + (1|gene), data = data)
#' base_model <- lmer(y ~ x + (1|gene), data = data)
#' linear_model <- lmer(y ~ poly(x, 2) + (1|gene), data = data)
#' cubic_model <- lmer(y ~ poly(x, 3) + (1|gene), data = data)
#' models_list <- list(null = null_model, base = base_model, linear = linear_model, cubic = cubic_model)
#' best <- select_best_model(models_list)
#' print(best)
#' }
select_best_model <- function(models) {
  if (!all(c("null","linear","cubic") %in% names(models))) {
    stop("The models list should contain 'null', 'linear', and 'cubic' models.")
  }

  # Compute AICc values
  aicc_values <- sapply(models, AICc)

  # Determine the best model based on criteria
  valid_models <- aicc_values[aicc_values <= (aicc_values["null"] - 2)]
  best_model <- names(which.min(valid_models))

  if (is.null(best_model)) {
    best_model <- "null"
  }

  return(best_model)
}
