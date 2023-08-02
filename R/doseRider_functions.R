#############
## Create SE Protein
###########

#library(gratia)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#########
## Genesets
###########
find_geneset_index <- function(gmt, geneset_name) {
  for (i in seq_along(gmt)) {
    if (gmt[[i]]$pathway == geneset_name) {
      return(i)
    }
  }
  return(NULL)
}

filter_gmt_by_size <- function(gmt, minGenesetSize, maxGenesetSize) {
  # Create an empty list to store filtered gene sets
  filtered_gmt <- list()

  # Loop over all gene sets in the GMT file
  for (i in seq_along(gmt)) {
    geneset <- gmt[[i]]$genes

    # Check if the number of genes in the gene set is within the specified range
    if (length(geneset) >= minGenesetSize & length(geneset) <= maxGenesetSize) {
      filtered_gmt[[i]] <- gmt[[i]]
    }
  }

  # Remove NULL elements from the list (i.e., the gene sets that were not within the specified range)
  filtered_gmt <- filtered_gmt[!sapply(filtered_gmt, is.null)]

  return(filtered_gmt)
}


load_consensupathdb_genesets <- function(filepath) {
  # Read the file and split into rows
  gene_sets_file <- readLines(filepath, warn = FALSE)

  # Process the gene sets file
  gene_sets <- gene_sets_file %>%
    strsplit("\t") %>%
    purrr::map(function(x) {
      list(
        pathway = x[1],
        external_id = x[2],
        source = x[3],
        genes = unlist(strsplit(x[4:length(x)], ","))
      )
    })

  return(gene_sets)
}

# Function to compute BgRatio and GeneRatio
## https://www.biostars.org/p/220465/
compute_ratios <- function(geneset, all_genes, input_vector) {
  M <- length(geneset)
  N <- length(all_genes)
  bgRatio <- M / N

  overlap_geneset <- length(intersect(input_vector, geneset))
  overlap_all <- length(intersect(input_vector, all_genes))
  geneRatio <- overlap_geneset / overlap_all

  return(list(bgRatio = bgRatio, geneRatio = geneRatio))
}


############
## Prepare data functions
#########
create_summarized_experiment <- function(assay_df, metadata_df) {
  # Convert metadata and assay data frames to appropriate objects
  coldata <- as.data.frame(metadata_df)
  assay <- as.matrix(assay_df)

  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = assay), colData = coldata)

  return(se)
}



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


################
## RNASeq Parameters
####################

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



###########
## Modelling
##########


# Function to fit GAM and compute AIC and BIC
fit_gam <- function(formula, data) {
  #data$start_event <- data$dose == 0
  tryCatch(
    {
      gam_model <- bam(as.formula(formula), data = data, method = "ML")
      return(gam_model)
    },
    error = function(e) {
      return(NA)
    }
  )
}



create_gamm_formula <- function(response, fixed_effects, random_effects, covariates = c(), model_type = "base") {
  format_random_effects_intercept <- paste0("s(", random_effects, ", bs = 're')")
  format_random_effects_slope <- paste0("s(",random_effects,",", fixed_effects, ", bs = 're')")
  format_fixed_effects <- paste0("s(", fixed_effects, ", bs = 'cr', k = 7)")

  base_formula <- paste(response, "~",
                        paste(format_random_effects_intercept, collapse = " + "),
                        if (length(covariates) > 0) paste(covariates, collapse = " + "))

  if (model_type == "base") {
    return(base_formula)
  } else if (model_type == "linear") {
    linear_formula <- paste(response, "~ ", paste(fixed_effects, collapse = " + "), "+",
                            paste(format_random_effects_intercept, collapse = " + "),
                            if (length(covariates) > 0) paste(covariates, collapse = " + "))
    return(linear_formula)
  } else if (model_type == "cubic") {
    cubic_formula <- paste(response, "~ ", paste(fixed_effects, collapse = " + "), "+",
                           paste(format_random_effects_intercept, collapse = " + "), "+", paste(format_fixed_effects, collapse = " + "),
                           "+", paste(format_random_effects_slope, collapse = " + ") ,
                           if (length(covariates) > 0) paste(covariates, collapse = " + "))
    return(cubic_formula)
  } else {
    stop("Invalid model type. Available options: 'base', 'linear', 'cubic'")
  }
}

compute_metrics <- function(model) {
  if (is.list(model)) {
    return(list(
      AIC = AIC(model),
      BIC = BIC(model),
      edf = sum(influence(model))
    ))
  } else {
    return(list(
      AIC = NA,
      BIC = NA,
      edf = NA
    ))
  }
}


################
## Main function
################

perform_analysis <- function(se, gmt, base_formula, linear_formula, cubic_formula, dose_col, sample_col, omic = "rnaseq") {
  # Create an empty data frame to store the results
  res_df <- data.frame(
    Geneset = character(),
    GeneRatio = numeric(),
    bgRatio = numeric(),
    Geneset_Size = integer(),
    Genes = integer(),
    Base_AIC = numeric(),
    Base_BIC = numeric(),
    Base_edf = numeric(),
    Linear_AIC = numeric(),
    Linear_BIC = numeric(),
    Linear_edf = numeric(),
    Cubic_AIC = numeric(),
    Cubic_BIC = numeric(),
    Cubic_edf = numeric(),
    P_Value_Linear = numeric(),
    P_Value_Cubic = numeric(),
    stringsAsFactors = FALSE
  )

  # Get the total number of gene sets
  total_gene_sets <- length(gmt)

  # Vector of uniques genes in the db
  all_genes <- unique(unlist(lapply(gmt, `[[`, "genes")))
  input_vector <- unique(rownames(se))
  # Initialize the progress bar
  pb <- progress::progress_bar$new(total = total_gene_sets)
  # Loop over the gene sets
  for (i in seq_along(gmt)) {
    geneset <- gmt[[i]]$genes

    # Filter gene set by length
    long_df <- prepare_data(se, geneset, dose_col, sample_col, omic)

    if (!is.null(long_df) && (length(long_df$gene) > 2)) {
      #Compute GeneRatio
      #input_vector <- unique(long_df$gene)
      ratios <- compute_ratios(geneset, all_genes, input_vector)

      # Fit GAM models and compute AIC and BIC
      # Fit GAM models and compute AIC, BIC, and df
      base_results <- fit_gam(base_formula, long_df)
      base_metrics <- compute_metrics(base_results)
      linear_results <- fit_gam(linear_formula, long_df)
      linear_metrics <- compute_metrics(linear_results)
      cubic_results <- fit_gam(cubic_formula, long_df)
      cubic_metrics <- compute_metrics(cubic_results)

      # Check if base_results and linear_results are not NaN before computing p-value
      if (is.list(base_results) && is.list(linear_results)) {
        p_value_linear <- as.numeric(itsadug::compareML(base_results, linear_results, signif.stars = FALSE, print.output = FALSE)$table[2, "p.value"])
      }

      # Check if base_results and cubic_results are not NaN before computing p-value
      if (is.list(base_results) && is.list(cubic_results)) {
        p_value_cubic <- as.numeric(itsadug::compareML(base_results, cubic_results, signif.stars = FALSE, print.output = FALSE)$table[2, "p.value"])
      }

      # Create a temporary data frame with the current results
      # Create a temporary data frame with the current results
      temp_df <- data.frame(
        Geneset = gmt[[i]]$pathway,
        Geneset_Size = length(geneset),
        GeneRatio = ratios$geneRatio,
        bgRatio = ratios$bgRatio,
        Genes = length(unique(long_df$gene)),
        Base_AIC = base_metrics$AIC,
        Base_BIC = base_metrics$BIC,
        Base_edf = base_metrics$edf,
        Linear_AIC = linear_metrics$AIC,
        Linear_BIC = linear_metrics$BIC,
        Linear_edf = linear_metrics$edf,
        Cubic_AIC = cubic_metrics$AIC,
        Cubic_BIC = cubic_metrics$BIC,
        Cubic_edf = cubic_metrics$edf,
        P_Value_Linear = p_value_linear,
        P_Value_Cubic = p_value_cubic,
        stringsAsFactors = FALSE
      )

      # Append the current results to the main results data frame
      res_df <- rbind(res_df, temp_df)
    }
    # Update the progress bar
    pb$tick()
  }

  # Add FDR adjustment to the results
  res_df$FDR <- p.adjust(res_df$P_Value_Cubic, method = "fdr")

  return(res_df)
}

#######
## Ploting
#######

library(ggplot2)

plot_smooth <- function(results, long_df, dose_col = "dose") {

  # Generate predictions with confidence intervals
  new_data <- expand.grid(
    dose = seq(min(long_df[dose_col]), max(long_df[dose_col]), length.out = 50),
    gene = unique(long_df$gene),
    sample = unique(long_df$sample)
  )

  colnames(new_data) <- c(dose_col,"gene","sample")
  predictions <- predict(results, newdata = new_data, se.fit = TRUE)
  predictions <- cbind(new_data, predictions)
  #trend_change_points <- compute_trend_change_points(results, long_df)

  # Center and scale the expression
  #predictions$fit <- (predictions$fit - mean(predictions$fit)) / sd(predictions$fit)
  #mean_data <- aggregate(fit ~ gene + dose, data = predictions, FUN = mean)

  # Plot the smooth curve with confidence interval region
  p <- ggplot(predictions, aes_string(x = dose_col, y = "fit", color = "gene")) +
    geom_line(aes(group = gene), color = "blue") +
    geom_smooth(formula = y ~ poly(x, 3), data = predictions, aes_string(x = dose_col, y = "fit"), se = TRUE, method = "lm", color = "red") +
    labs(x = "Dose", y = "Normalized Expression") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),     # Increase the text size of axis labels
          axis.title = element_text(size = 14),    # Increase the text size of axis titles
          plot.title = element_text(size = 16))


  return(p)
}



plot_raw <- function(results, long_df) {
  # Generate predictions with confidence intervals


  # Center and scale the expression
  #predictions$fit <- (predictions$fit - mean(predictions$fit)) / sd(predictions$fit)
  mean_data <- aggregate(counts ~ gene + dose, data = long_df, FUN = mean)

  # Plot the smooth curve with confidence interval region
  p <- ggplot(mean_data, aes(x = dose, y = counts, color = gene)) +
    geom_line(aes(group = gene), color = "blue") +
    geom_smooth(formula = y ~ poly(x, 3), data = mean_data, aes(x = dose, y = counts), se = F, method = "lm", color = "red") +
    labs(x = "Dose", y = "Normalized Expression") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12),     # Increase the text size of axis labels
          axis.title = element_text(size = 14),    # Increase the text size of axis titles
          plot.title = element_text(size = 16))


  return(p)
}

compute_trend_change_points <- function(gamm_model, long_df) {
  # Get the predicted values from the model
  fixed_smooth <- predict(gamm_model, type = "terms")[,1]

  # Calculate the dose at which the maximum change in gene expression occurs
  global_max_dose <- long_df$dose[which.max(fixed_smooth)]

  # Calculate the dose at which the minimum change in gene expression occurs
  global_min_dose <- long_df$dose[which.min(fixed_smooth)]

  # Calculate the first derivative of the predicted values
  derivative <- diff(fixed_smooth)

  # The points where the derivative changes sign are local maxima and minima
  local_max_min_indices <- which(diff(sign(derivative)) != 0)

  # Get the local maxima and minima doses
  local_max_min_doses <- long_df$dose[local_max_min_indices]

  # Create vectors for local max and min
  local_max_doses <- local_max_min_doses[which(fixed_smooth[local_max_min_indices] > 0)]
  local_min_doses <- local_max_min_doses[which(fixed_smooth[local_max_min_indices] < 0)]

  return(list(global_max_dose = unique(global_max_dose),
              global_min_dose = unique(global_min_dose),
              local_max_doses = unique(local_max_doses),
              local_min_doses = unique(local_min_doses)))
}

