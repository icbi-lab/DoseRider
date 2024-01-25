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
#' @importFrom mgcv bam
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
fit_gam <- function(formula, data, omic) {
  # Adjust control parameters
  #control_params <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e7))

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
        gam_model <- bam(as.formula(formula), data = data, family = family,  method = "ML")
      }
      return(gam_model)
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
create_gamm_formula <- function(response, fixed_effects, random_effects, covariates = c(), model_type = "non_linear", omic = NULL, k = NULL) {
  base_formula <- paste(response, "~")

  if (model_type == "null") {
    # Null model with only random intercept
    base_formula <- paste(base_formula, "s(", random_effects, ", bs='re')")
  } else if (model_type == "linear") {
    # Linear model with random slope and intercept
    base_formula <- paste(base_formula, fixed_effects, "+ s(", random_effects, ", bs='re')", "+ s(",fixed_effects, ",",random_effects,",bs='re')")
  } else if (model_type == "non_linear") {
    # Non-linear model with random intercept and smooth terms
    if (is.null(k)) k <- 3 # Default knots if not specified
    #base_formula <- paste(base_formula, " s(", fixed_effects, ", bs='cr', k=", k, ")", "+ s(", fixed_effects, ",", random_effects, ", bs='fs', k =",k,")")
    base_formula <- paste(base_formula, " s(", fixed_effects, ", bs='cr', k=", k, ")", "+ s(", random_effects, ", bs='re')", "+ s(",fixed_effects, ",",random_effects,",bs='re')")

  }

  # Add covariates if any
  if (length(covariates) > 0) {
    covariate_terms <- paste(covariates, collapse = " + ")
    base_formula <- paste(base_formula, "+", covariate_terms)
  }

  # Add offset for RNA-Seq data
  if (!is.null(omic) && omic == "rnaseq") {
    base_formula <- paste(base_formula, "+ offset(log(size_factor))")
  }

  return(base_formula)
}



#' Extract Random Effects from GAMM
#'
#' This function extracts the random effect estimates for each gene from a GAMM object.
#' The result is a dataframe with genes as rows and their corresponding random effect values.
#'
#' @param model A GAMM model object.
#' @return A dataframe with genes as rows and their random effect estimates.
#' @importFrom gammit extract_ranef
extract_random_effects_gamm <- function(model, dose_col) {
  # Ensure the input is a GAMM model
  if (!inherits(model, c("gam", "bam"))) {
    stop("The provided model is not a GAMM model.")
  }

  # Extract the coefficients from the model
  random_effects <- extract_ranef(model)
  random_slope <- as.data.frame(random_effects[random_effects$effect == dose_col,])
  random_intercept <- as.data.frame(random_effects[random_effects$effect == "Intercept",])
  rownames(random_slope) <- random_slope$group
  rownames(random_intercept) <- random_intercept$group

  # Convert the random effects into a dataframe
  random_effects_df <- data.frame(row.names = rownames(random_slope), RandomIntercept = unlist(random_intercept$value, use.names = FALSE), RandomEffect = unlist(random_slope$value, use.names = FALSE))

  return(random_effects_df)
}



#' Compute metrics for a given gam or bam model
#'
#' This function computes the AIC, BIC, and effective degrees of freedom (edf) for a given `gam` or `bam` model.
#'
#' The approach to calculate the effective degrees of freedom (edf)
#' is adapted from the itsadug package: Interpreting Time Series and Autocorrelated Data Using GAMMs
#' Reference: itsadug package (https://CRAN.R-project.org/package=itsadug)
#'
#' @param model A `gam` or `bam` model object from the mgcv package.
#'
#' @return A list containing the AIC, BIC, and edf of the model.
#'         If the model parameter isn't a `gam` or `bam` object, the function will return NA for all metrics.
#' @examples
#' library(mgcv)
#' data(sleepstudy)
#' gam_model <- gam(Reaction ~ s(Days), data = sleepstudy)
#'
#' # Compute metrics
#' compute_metrics_gamm(gam_model)
#' @importFrom AICcmodavg AICc
compute_metrics_gamm <- function(model) {
  if (inherits(model, c("gam", "bam"))) {
    return(list(
      AIC = AIC(model),
      AICc = AICc(model),
      BIC = BIC(model),
      edf = length(model$sp) + model$nsdf + ifelse(length(model$smooth) > 0,
            sum(sapply(model$smooth, FUN = function(x) {x$null.space.dim},
                                                        USE.NAMES = FALSE)), 0)
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
