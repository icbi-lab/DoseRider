#' Create a SummarizedExperiment object from assay and metadata data frames
#'
#' This function creates a SummarizedExperiment object which combines the assay
#' data (counts, or other measurements) and associated metadata.
#'
#' @param assay_df Data frame containing assay data, with rows as genes and columns as samples.
#' @param metadata_df Data frame containing metadata for the samples.
#'
#' @return A SummarizedExperiment object with the given assay and metadata.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @examples
#' assay_df <- data.frame(matrix(runif(100), nrow=10))
#' metadata_df <- data.frame(condition = rep(c("A", "B"), each = 5))
#' se <- create_summarized_experiment(assay_df, metadata_df)
#' @export
create_summarized_experiment <- function(assay_df, metadata_df) {
  # Convert metadata and assay data frames to appropriate objects
  coldata <- as.data.frame(metadata_df)
  assay <- as.matrix(assay_df)

  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = assay), colData = coldata)

  return(se)
}

#' Prepare data for use with glmmTMB from a SummarizedExperiment object
#'
#' This function prepares the data from a SummarizedExperiment object into a format
#' suitable for the glmmTMB function. This includes transforming the data to long format,
#' adding size factors and dispersions (if applicable), and merging with sample metadata.
#'
#' @param se A SummarizedExperiment object containing both assay and metadata.
#' @param geneset A character vector of gene names to be considered.
#' @param dose_col The name of the dose column in the metadata of the SummarizedExperiment object.
#' @param sample_col The name of the sample column in the metadata of the SummarizedExperiment object.
#' @param omic A character string specifying the type of omics data. Default is "rnaseq".
#'
#' @return A data frame in a format suitable for glmmTMB function.
#' @importFrom reshape melt
#'
#' @examples
#' geneset <- rownames(assay_df)
#' dose_col <- "condition"
#' sample_col <- "sample"
#' long_df <- prepare_data(se, geneset, dose_col, sample_col)
#' @export

prepare_data <- function(se, geneset, dose_col, sample_col, omic = "rnaseq") {

  # Obtain size_factors
  size_factors <- colData(se)$size_factors

  # Prepare data for glmmTMB, including size factors and dispersions
  common_genes <- intersect(geneset, rownames(se))
  if (length(common_genes) > 0) {
    long_df <- suppressWarnings(as.data.frame(reshape::melt(t(assay(se)[common_genes,]), as.is = TRUE), warning = FALSE))
    colnames(long_df) <- c("sample", "gene", "counts")
    long_df$size_factor <- size_factors[match(long_df$sample, rownames(colData(se)))]
    long_df$dispersion <- colData(se)$dispersion

    if (omic != "rnaseq") {
      long_df$size_factor <- NA
      long_df$dispersion <- NA
    }

    long_df <- merge(long_df, colData(se), by = "sample")
    long_df$dose <- unlist(as.numeric(long_df[[dose_col]]))
    long_df <- as.data.frame(long_df, warning = FALSE)
  } else {
    long_df <- NULL
  }

  return(long_df)
}

#' Estimate model parameters for DESeq2 or edgeR
#'
#' This function estimates necessary parameters (size factors and dispersions) for either DESeq2 or edgeR
#' based on the input count data and sample metadata.
#'
#' @param count_data A matrix or data frame of count data, with rows as genes and columns as samples.
#' @param sample_metadata A data frame of sample metadata matching the columns of count_data.
#' @param model_type A character string specifying the type of model. Must be either "DESeq2" or "edgeR".
#'
#' @return A list containing the estimated size factors, dispersions, and theta (1/dispersion).
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions
#' @importFrom edgeR DGEList calcNormFactors estimateDisp
#'
#' @examples
#' sample_metadata <- data.frame(condition = rep(c("A", "B"), each = 5))
#' parameters <- estimate_model_parameters(assay_df, sample_metadata, model_type = "DESeq2")

estimate_model_parameters <- function(count_data, sample_metadata, model_type) {
  if (model_type == "DESeq2") {
    library(DESeq2)

    # Create DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(countData = count_data,
                                  colData = sample_metadata,
                                  design = ~ dose)

    # Estimate size factors
    dds <- estimateSizeFactors(dds)
    size_factors <- sizeFactors(dds)

    # Estimate dispersions
    dds <- estimateDispersions(dds)
    dispersions <- dispersions(dds)

    # Set theta to 1/dispersions
    theta <- 1 / dispersions

  } else if (model_type == "edgeR") {
    library(edgeR)

    # Create DGEList object
    dge <- DGEList(counts = count_data, group = sample_metadata$dose)

    # Estimate normalization factors
    dge <- calcNormFactors(dge)
    size_factors <- dge$samples$norm.factors

    # Estimate dispersion
    dge <- estimateDisp(dge)
    dispersions <- dge$common.dispersion

    # Set theta to 1/dispersions
    theta <- 1 / dispersions


  } else {
    stop("Invalid model_type specified. Must be either 'DESeq2' or 'edgeR'.")
  }

  # Return a list containing the estimated parameters
  parameters <- list(dispersions = dispersions,
                     size_factors = size_factors,
                     theta = theta)

  return(parameters)
}

