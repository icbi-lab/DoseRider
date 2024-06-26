% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_models.R
\name{compare_all_models}
\alias{compare_all_models}
\title{Compare Models with ANOVA}
\usage{
compare_all_models(
  null_results,
  linear_results,
  non_linear_fixed_results,
  non_linear_mixed_results,
  modelType
)
}
\arguments{
\item{linear_results}{A fitted model representing the linear model.}

\item{base_results}{A fitted model representing the base model.}

\item{cubic_results}{A fitted model representing the cubic model.}
}
\value{
A named vector of p-values for the comparisons: base vs. linear and linear vs. cubic.
}
\description{
This function compares the goodness of fit of three models (base, linear, and cubic)
using an ANOVA test and returns the p-values associated with the model comparisons.
}
\examples{
# Assuming lm_base, lm_linear, and lm_cubic are fitted linear models:
# compare_models_pvalues <- compare_three_models(lm_base, lm_linear, lm_cubic)
}
