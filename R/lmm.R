#' Fit Linear Mixed Model (LMM) or Generalized Linear Mixed Model (GLMM)
#'
#' This function fits an LMM or a GLMM to the provided data using the specified formula.
#' The aim is to model all the genes from a specified pathway to check the effects of a dose.
#' The model's family is chosen based on the 'omic' parameter.
#'
#' @param formula A formula for the LMM or GLMM.
#' @param data A data frame (in long format) containing the data to be modeled. The data frame should include
#' columns for dose, gene (or metabolite), and, if 'rnaseq' is specified, an offset. Other covariates may also be present.
#' @param omic A character string specifying the type of data. If 'rnaseq', the family is set to negative binomial.
#' @param theta The theta parameter for the negative binomial family. Relevant only if omic = "rnaseq".
#'
#' @return An LMM or GLMM model if fitting is successful; NA otherwise.
#'
#' @importFrom lme4 glmer lmer glmerControl
#' @importFrom MASS negative.binomial
#' @importFrom stats formula
#' @importFrom utils globalVariables
#' @importFrom splines bs
#'
#' @examples
#' \dontrun{
#' data("mtcars")
#' # Generate a suitable formula for the LMM
#' formula <- mpg ~ hp + (1|cyl)
#' model <- fit_lmm(formula = formula, data = mtcars, omic = "rnaseq", theta = 1.5)
#' }
#'
fit_lmm <- function(formula, data, omic = "base") {
  # Adjust control parameters
  control_params <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e8))

  # Determine the family based on 'omic' parameter
  if (omic == "rnaseq") {
    #if (!"theta" %in% colnames(data)) stop("theta must be specified for omic = 'rnaseq'")
    family <- negative.binomial(theta = unique(data$theta))
    family_choice <- "negative_binomial"
  } else {
    family_choice <- "gaussian"
  }

  tryCatch(
    {
      if (family_choice == "gaussian") {
        lmm_model <- lmer(as.formula(formula), data = data, REML = F)
      } else {
        lmm_model <- glmer(as.formula(formula), data = data, family = family, control = control_params)
      }
      return(lmm_model)
    },
    error = function(e) {
      cat(paste0("\n[!] Error: ",e))
      return(NA)
    }
  )
}



#' Create Generalized Linear Mixed Model (GLMM) Formula
#'
#' This function generates a formula for a GAMM, given response and effect variables.
#' The main goal is to compare the 'cubic' model with the 'base' model.
#'
#' @param response A string specifying the response variable in the model.
#' @param fixed_effects A string or a vector of strings specifying the fixed effects in the model.
#' @param random_effects A string or a vector of strings specifying the random effects in the model.
#' @param covariates A string or a vector of strings specifying any covariates to be included in the model. Default is an empty vector.
#' @param model_type A string specifying the type of the model. Can be 'base' or 'cubic'. Default is 'base'.
#' @param omic A character string specifying the type of data. If 'rnaseq', an offset is included in the formula.
#'
#' @return A string representing the formula for the GAMM.
#'
#' @examples
#' \dontrun{
#' # Create a base formula
#' base_formula <- create_gamm_formula(response = "counts",
#'                                     fixed_effects = "dose",
#'                                     random_effects = "gene",
#'                                     model_type = "base",
#'                                     omic = "rnaseq")
# # Create a cubic formula
#' cubic_formula <- create_gamm_formula(response = "counts",
#'                                      fixed_effects = "dose",
#'                                      random_effects = "gene",
#'                                      model_type = "cubic",
#'                                      omic = "rnaseq")
#' }
create_lmm_formula <- function(response, fixed_effects, random_effects, covariates = c(), model_type = "base", omic = NULL, k = NULL) {

  base_formula <- paste(response, "~")

  if (model_type == "null") {
    base_formula <- paste(base_formula, "(1|", random_effects, ")")
    if (length(covariates) > 0) {
      base_formula <- paste(base_formula, "+", paste(covariates, collapse = " + "))
    }
    if (!is.null(omic) && omic == "rnaseq") {
      base_formula <- paste(base_formula, "+ offset(log(size_factor))")
    }
  } else {
    base_formula <- paste(base_formula, "(", fixed_effects, "|", random_effects, ")")

    if (!is.null(omic) && omic == "rnaseq") {
      base_formula <- paste(base_formula, "+ offset(log(size_factor))")
    }

    if (length(covariates) > 0) {
      base_formula <- paste(base_formula, "+", paste(covariates, collapse = " + "))
    }

    if (model_type == "cubic") {
      cubic_terms <- paste("bs(", fixed_effects, ",k=", k -1, ")", collapse = " + ")
      base_formula <- paste(base_formula, "+", cubic_terms)
    } else if (model_type == "linear") {
      linear_terms <- fixed_effects
      base_formula <- paste(base_formula, "+", linear_terms)
    } else if (model_type != "base") {
      stop("Invalid model type. Available options: 'base', 'linear', 'cubic', 'null'")
    }
  }

  return(base_formula)
}

#' Extract Random Effects from LMM
#'
#' This function extracts the random effect estimates for each gene from an LMM model object.
#' The result is a dataframe with genes as rows and their corresponding random effect values.
#'
#' @param model An LMM model object, typically of class 'lmerMod'.
#' @return A dataframe with genes as rows and their random effect estimates.
#' @importFrom lme4 ranef
extract_random_effects_lmm <- function(model, dose_col) {
  # Ensure the input is an LMM model
  if (!inherits(model, "merMod")) {
    stop("The provided model is not an LMM model.")
  }

  # Extract the random effects using ranef
  random_effects <- ranef(model)$gene

  # If there are multiple random effects (e.g., intercept and slope), you may need to adjust this part
  # Here, we assume a single random effect for simplicity
  random_effects_df <- data.frame(row.names = rownames(random_effects), RandomEffect = random_effects[,dose_col])

  return(random_effects_df)
}


#' Compute metrics for a given lmer or glmer model
#'
#' This function computes the AIC, AICc, BIC, and degrees of freedom metrics for a given `lmer` or `glmer` model.
#'
#' @param model A `lmer` or `glmer` model object. If it's neither, the function returns NA for all the metrics.
#'
#' @return A list containing the AIC, BIC, and degrees of freedom of the model. If the model parameter isn't a `lmer` or `glmer` object, the function will return NA for all metrics.
#' @examples
#' library(lme4)
#' data(sleepstudy)
#' model <- lmer(Reaction ~ Days + (1|Subject), data = sleepstudy)
#'
#' # Compute metrics
#' compute_metrics(model)
#' @importFrom AICcmodavg AICc
compute_metrics_lmm <- function(model) {
  if (inherits(model, "merMod")) { # Checking if the model is either a lmer or glmer
    df_model <- df.residual(model) # degrees of freedom for the model
    return(list(
      AIC = AIC(model),
      AICc = AICc(model),
      BIC = BIC(model),
      edf = df_model
    ))
  } else {
    return(list(
      AIC = NA,
      AICc = NA,
      BIC = NA,
      edf = NA
    ))
  }
}


