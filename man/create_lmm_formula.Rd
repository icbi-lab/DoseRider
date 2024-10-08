% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmm.R
\name{create_lmm_formula}
\alias{create_lmm_formula}
\title{Create Generalized Linear Mixed Model (GLMM) Formula}
\usage{
create_lmm_formula(
  response,
  fixed_effects,
  random_effects,
  covariates = c(),
  model_type = "base",
  omic = NULL,
  spline_knots = 3
)
}
\arguments{
\item{response}{A string specifying the response variable in the model.}

\item{fixed_effects}{A string or a vector of strings specifying the fixed effects in the model.}

\item{random_effects}{A string or a vector of strings specifying the random effects in the model.}

\item{covariates}{A string or a vector of strings specifying any covariates to be included in the model. Default is an empty vector.}

\item{model_type}{A string specifying the type of the model. Can be 'base' or 'cubic'. Default is 'base'.}

\item{omic}{A character string specifying the type of data. If 'rnaseq', an offset is included in the formula.}
}
\value{
A string representing the formula for the GAMM.
}
\description{
This function generates a formula for a GAMM, given response and effect variables.
The main goal is to compare the 'cubic' model with the 'base' model.
}
\examples{
\dontrun{
# Create a base formula
base_formula <- create_gamm_formula(response = "counts",
                                    fixed_effects = "dose",
                                    random_effects = "gene",
                                    model_type = "base",
                                    omic = "rnaseq")
cubic_formula <- create_gamm_formula(response = "counts",
                                     fixed_effects = "dose",
                                     random_effects = "gene",
                                     model_type = "cubic",
                                     omic = "rnaseq")
}
}
