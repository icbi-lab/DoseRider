#
# library(dplyr)
# load_consensupathdb_genesets <- function(filepath) {
# # Read the file and split into rows
# gene_sets_file <- readLines(filepath, warn = FALSE)
#
# # Process the gene sets file
# gene_sets <- gene_sets_file %>%
#   strsplit("\t") %>%
#   purrr::map(function(x) {
#     list(
#       pathway = x[1],
#       external_id = x[2],
#       source = x[3],
#       genes = unlist(strsplit(x[4:length(x)], ","))
#     )
#   })
#
# return(gene_sets)
# }
#
# # # List of all file names
# # file_names <- c("CPDB_pathways_genes_ensembl.tab",
# #                 "CPDB_pathways_genes_entrez.tab", "CPDB_pathways_genes_refseq.tab",
# #                 "CPDB_pathways_metabolites.tab")
#
#
# CPDB_pathways_genes_ensembl <- load_consensupathdb_genesets("../data/CPDB_pathways_genes_ensembl.tab")
# CPDB_pathways_genes_entrez<- load_consensupathdb_genesets("../data/CPDB_pathways_genes_entrez.tab")
# CPDB_pathways_genes_refseq <- load_consensupathdb_genesets("../data/CPDB_pathways_genes_refseq.tab")
# CPDB_pathways_genes_symbol <- load_consensupathdb_genesets("../data/CPDB_pathways_genes_symbol.tab")
# CPDB_pathways_metabolites <- load_consensupathdb_genesets("../data/CPDB_pathways_metabolites.tab")
#
# CPDB_pathways_genes_ensembl <- CPDB_pathways_genes_ensembl
#
# usethis::use_data(CPDB_pathways_genes_ensembl, overwrite = T)
# usethis::use_data(CPDB_pathways_genes_entrez, overwrite = T)
# usethis::use_data(CPDB_pathways_genes_symbol, overwrite = T)
# usethis::use_data(CPDB_pathways_genes_refseq,overwrite = T)
# usethis::use_data(CPDB_pathways_metabolites, overwrite = T)
