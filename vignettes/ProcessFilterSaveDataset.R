library(doseRider)
library(DESeq2)
library(ggplot2)

# Load the dataset
load("data/PRJNA869442.rda")

# Define the compound doses
compound_doses <- list(
  "2,4-BPF" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "2,4-BPS" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
  "4,4-BPF" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50),
  "BADGE" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "BPA" = c(0.0005, 0.1, 0.5, 1, 5, 10),
  "BPAF" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "BPAP" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5),
  "BPC" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "BPS" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
  "BPS-MAE" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "BPS-MPE" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "BTUM" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "D-8" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10),
  "DDS" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5),
  "P201" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50),
  "TGSA" = c(0.0005, 0.001, 0.01, 0.1, 0.5, 1, 5, 10)
)

# Function to process each compound
process_compound <- function(compound, doses) {
  print(paste0("[+] Analyzing Compound: ",compound))
  # Filter data for the current compound and the selected doses
  # Example filtering for BPA
  compound_data <- PRJNA869442[,(colData(PRJNA869442)$Chemical %in% c("Cells, no treatment",compound)) &
                             (colData(PRJNA869442)$Dose %in% doses)]

  compound_data$sample <- colnames(compound_data)

  # Differential expression analysis
  colData(compound_data)$Dose <- as.numeric(colData(compound_data)$Dose)

  dds <- DESeqDataSetFromMatrix(countData = assay(compound_data),
                                colData = colData(compound_data),
                                design = ~ Dose * Dose)

  dds <- DESeq(dds, parallel = TRUE, quiet = TRUE)
  res <- results(dds)
  res_df <- as.data.frame(res)

  # Filter significant genes
  filter_res <- res_df[(res_df$baseMean > 11) & (res_df$pvalue < 0.1), ]
  compound_data <- estimate_model_parameters(compound_data)

  filter_compound_data <- compound_data[rownames(compound_data) %in% rownames(filter_res), ]

  # Extract metadata and expression data
  metadata <- colData(filter_compound_data)
  expression_data <- assay(filter_compound_data)
  if (!dir.exists("../PRJNA869442/")) {
    dir.create("../PRJNA869442/",showWarnings = F,recursive = T)

  }
  # Save metadata and expression data to TSV files
  metadata_file <- paste0("../PRJNA869442/",compound, "_metadata.tsv")
  expression_file <- paste0("../PRJNA869442/",compound, "_expression_data.tsv")

  write.table(metadata, metadata_file, sep = "\t", quote = FALSE)
  write.table(expression_data, expression_file, sep = "\t", quote = FALSE)
}

# Loop through each compound and process it
for (compound in names(compound_doses)) {
  process_compound(compound, c(0, compound_doses[[compound]]))
}

