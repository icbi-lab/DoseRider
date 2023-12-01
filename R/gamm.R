#' Fit Generalized Additive Model (GAM)
#'
#' This function fits a GAM to the provided data using the specified formula.
#' It adapts to different omics types and handles family specification for RNAseq data.
#'
#' @param formula A formula for the GAM.
#' @param data A data frame containing the data for modeling.
#' @param omic A character string specifying the type of omic data.
#'             If 'rnaseq', a negative binomial family is used.
#'
#' @return A GAM model if fitting is successful; NA otherwise.
#'
#' @importFrom mgcv gam
#' @importFrom stats formula
#' @importFrom MASS negative.binomial
#' @importFrom utils globalVariables
#'
#' @examples
#' \dontrun{
#'   data("mtcars")
#'   formula <- mpg ~ s(hp) + s(wt)
#'   model <- fit_gam(formula, data = mtcars, omic = "base")
#' }
#'
#' @export
fit_gam <- function(formula, data) {
  # Adjust control parameters
  control_params <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e7))
  
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
        gam_model <- bam(as.formula(formula), data = data, method = "ML")
      } else {
        gam_model <- bam(as.formula(formula), data = data, family = family, control = control_params,  method = "ML")
      }
      return(lmm_model)
    },
    error = function(e) {
      cat(paste0("[!] Error: ",e))
      return(NA)
    }
  )
}


#' Create Generalized Additive Mixed Model (GAMM) Formula with Omics Consideration
#'
#' Generates a formula for fitting a GAMM, given response, fixed effects, random effects, omics type, and model type.
#' Supports 'base', 'linear', and 'cubic' models and can handle RNA-seq data.
#'
#' @param response A string specifying the response variable in the model.
#' @param fixed_effects A string or a vector of strings specifying the fixed effects in the model.
#' @param random_effects A string or a vector of strings specifying the random effects in the model.
#' @param covariates Optional vector of additional covariates in the model.
#' @param model_type Type of the model ('base', 'linear', 'cubic').
#' @param omic Type of omic data ('base' or 'rnaseq').
#'
#' @return A string representing the GAMM formula.
#'
#' @examples
#' \dontrun{
#'   base_formula <- create_gamm_formula("mpg", "hp", "cyl", model_type = "base", omic = "base")
#'   cubic_formula <- create_gamm_formula("mpg", "hp", "cyl", model_type = "cubic", omic = "rnaseq")
#' }
#' @export
create_gamm_formula <- function(response, fixed_effects, random_effects, covariates = c(), model_type = "base", omic = "base") {
  # Basic structure for random effects
  random_effect_structure <- paste0("s(", random_effects, ", bs = 're')")
  
  # Start assembling the base formula
  base_formula <- paste0(response, " ~ ", random_effect_structure)
  
  # Handling covariates
  if (length(covariates) > 0) {
    covariate_terms <- paste(covariates, collapse = " + ")
    base_formula <- paste(base_formula, "+", covariate_terms)
  }
  
  # Add omic-specific considerations
  if (omic == "rnaseq") {
    base_formula <- paste(base_formula, "+ offset(log(size_factor))")
  }
  
  # Extend formula based on model type
  if (model_type == "linear") {
    linear_terms <- paste(fixed_effects, collapse = " + ")
    return(paste(base_formula, "+", linear_terms))
  } else if (model_type == "cubic") {
    cubic_terms <- paste0("s(", fixed_effects, ", bs = 'cr', k = 4)") # Adjust 'k' as needed
    return(paste(base_formula, "+", cubic_terms))
  } else if (model_type != "base") {
    stop("Invalid model type. Available options: 'base', 'linear', 'cubic'")
  }
  
  return(base_formula)
}



