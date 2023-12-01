
#' Compare Two GAM Models
#'
#' This function compares two Generalized Additive Models (GAMs) by performing a Chi-square test on their scores.
#' The function first checks if both models are GAMs and if they are built using the same method.
#' It then computes the test statistic based on the difference of the models' Un-Biased Risk Estimators (UBREs)
#' and their degrees of freedom.
#'
#' @param model1 A gam or bam model object. The first model to compare.
#' @param model2 A gam or bam model object. The second model to compare.
#'
#' @return The p-value from the Chi-square test comparing the two models.
#'
#' @examples
#' data(mtcars)
#' fit1 <- mgcv::gam(mpg ~ s(hp), data = mtcars)
#' fit2 <- mgcv::gam(mpg ~ s(hp, k = 5), data = mtcars)
#' compareGAMM(fit1, fit2)
#'
#' @importFrom mgcv gam
#' @importFrom stats pchisq
#' @export
compareGAMM <- function(model1, model2) {
  # Check if models are gam or bam objects:
  if (!("gam" %in% class(model1)) || !("gam" %in% class(model2))) {
    stop("Models are not gam objects (i.e., built with bam()/gam()).")
  }

  # Check if models are comparable:
  if (model1$method != model2$method) {
    stop(sprintf("Models are incomparable: method model1 = %s, method model2 = %s", model1$method, model2$method))
  }

  # Extract the number of degrees of freedom (ndf) for each model:
  ndf1 <- length(model1$sp) + model1$nsdf + ifelse(length(model1$smooth) > 0,
                                                   sum(sapply(model1$smooth, FUN = function(x) {
                                                     x$null.space.dim
                                                   }, USE.NAMES = FALSE)), 0)
  ndf2 <- length(model2$sp) + model2$nsdf + ifelse(length(model2$smooth) > 0,
                                                   sum(sapply(model2$smooth, FUN = function(x) {
                                                     x$null.space.dim
                                                   }, USE.NAMES = FALSE)), 0)

  # Calculate AIC scores for the models:
  ml1 <- AIC(model1)
  ml2 <- AIC(model2)

  # Perform the chi-square test to obtain the p-value for model comparison:
  p_value <- pchisq(2 * (ml1 - ml2), abs(ndf1 - ndf2), lower.tail = FALSE)

  return(p_value)
}

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
#' @export
compare_all_models <- function(null_results, linear_results, cubic_results) {
  # Perform an ANOVA comparison between the three models
  anova_results <- anova(null_results, linear_results, cubic_results)

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
#' @export
select_best_model <- function(models) {
  if (!all(c("null","cubic") %in% names(models))) {
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
