% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doseRider.R
\name{DoseRider}
\alias{DoseRider}
\title{Perform DoseRider Analysis Using Linear Mixed Models (LMMs) or Generalized Additive Mixed Models (GAMMs)}
\usage{
DoseRider(
  se,
  gmt,
  dose_col = "dose",
  sample_col = "sample",
  covariates = c(),
  omic = "rnaseq",
  minGSsize = 5,
  maxGSsize = 300,
  method = "fdr",
  modelType = "LMM",
  FilterPathway = FALSE,
  pca_threshold = 0.6,
  log_transform = F,
  models = c("linear", "non_linear_fixed", "non_linear_mixed")
)
}
\arguments{
\item{se}{SummarizedExperiment object or a matrix/data frame containing gene expression data.}

\item{gmt}{List of gene sets, each represented as a list with gene names.}

\item{dose_col}{Name of the column representing dose information.}

\item{sample_col}{Name of the column representing sample information.}

\item{covariates}{Optional, vector specifying the covariate column(s) in `se`.}

\item{omic}{Type of omics data, defaults to "rnaseq".}

\item{minGSsize}{Minimum gene set size for analysis, defaults to 5.}

\item{maxGSsize}{Maximum gene set size for analysis, defaults to 300.}

\item{method}{Method for multiple testing adjustment, defaults to "fdr".}

\item{modelType}{Type of model to be used for analysis, "LMM" for Linear Mixed Models or "GAMM"
for Generalized Additive Mixed Models. Defaults to "LMM".}

\item{FilterPathway}{Boolean, if TRUE the function will apply PCA filtering to detect antagonist patterns. Defaults to FALSE.}

\item{pca_threshold}{Numeric value specifying the variance threshold for PC1 to filter pathways. Default is 0.6.}

\item{log_transform}{Logical, whether to log10 transform the dose values. Default is FALSE.}
}
\value{
A list containing results for each gene set including various metrics, p-values,
and adjusted p-values.
}
\description{
This function performs DoseRider analysis on gene expression data, applying either Linear Mixed Models (LMMs)
or Generalized Additive Mixed Models (GAMMs) to each gene set defined in the gene matrix transposed (GMT) format.
It evaluates dose-response relationships in the context of gene sets and calculates various model metrics,
significance, and smoothing predictions.
}
\examples{
\dontrun{
data("SummarizedExperiment")
gmt <- list(geneSet1 = list(genes = c("gene1", "gene2")))
results <- DoseRider(se, gmt, "dose", "sample", "covariate", "rnaseq", modelType = "GAMM")
}

}
