#' Find the Index of a Geneset in a GMT Object
#'
#' This function searches for a specified geneset within a GMT (Gene Matrix Transposed file format) object and returns its index.
#' If the geneset is not found, NULL is returned.
#'
#' @param gmt A list, where each element represents a geneset. Each geneset should be a named list with at least a field corresponding to the geneset name.
#' @param geneset_name A character string representing the name of the geneset to find.
#' @param pathway_col A character string representing the name of the field in the geneset list items which corresponds to the geneset name.
#'
#' @return The index of the geneset in the GMT object if found, otherwise NULL.
#' @examples
#' # Consider a GMT object
#' gmt <- list(list(pathway = "geneset1", gene_ids = c("gene1", "gene2")),
#'              list(pathway = "geneset2", gene_ids = c("gene3", "gene4")))
#'
#' # Find the index of a geneset
#' find_geneset_index(gmt, "geneset2", "pathway") # returns 2
#' find_geneset_index(gmt, "nonexistent_geneset", "pathway") # returns NULL
#'
#' @export
find_geneset_index <- function(gmt, geneset_name, pathway_col) {
  for (i in seq_along(gmt)) {
    if (gmt[[i]][[pathway_col]] == geneset_name) {
      return(i)
    }
  }
  return(NULL)
}

#' Read GMT File from MSigDB
#'
#' Parses a GMT (Gene Matrix Transposed) file from the Molecular Signatures Database (MSigDB),
#' structuring it into a list of pathways and their associated genes. Each element of the list
#' is a list itself, containing the pathway name and a vector of gene symbols.
#'
#' @param file_path The path to the GMT file to be read.
#' @return A list where each element is a list with two elements: `$pathway`, the name of the pathway,
#' and `$genes`, a vector of gene symbols associated with that pathway.
#' @examples
#' file_path <- "path/to/your/c2.cgp.v2023.2.Hs.symbols.gmt"
#' gmt_data <- read_gmt(file_path)
#' if (length(gmt_data) > 0) {
#'   print(gmt_data[[1]])
#'   print(gmt_data[[2]])
#' }
#' @export
read_gmt <- function(file_path) {
  con <- file(file_path, "r")
  gmt <- list()
  while(TRUE) {
    line <- readLines(con, n = 1, warn = FALSE)
    if(length(line) == 0) {
      break
    }
    elements <- strsplit(line, "\t")[[1]]
    pathway_name <- elements[1]
    genes <- elements[-c(1, 2)] # Exclude pathway name and URL
    gmt[[length(gmt) + 1]] <- list(pathway = pathway_name, genes = genes)
  }
  close(con)
  return(gmt)
}


#' Filter GMT by Geneset Size
#'
#' This function filters a GMT (Gene Matrix Transposed file format) object by a specified geneset size range.
#' Each geneset in the GMT object is checked to see if its size falls within the minimum and maximum size thresholds.
#' If it does, the geneset is kept; if it doesn't, the geneset is discarded.
#'
#' @param gmt A list, where each element represents a geneset. Each geneset should be a named list with a field for the geneset name and a field for the genes contained in the geneset.
#' @param minGenesetSize An integer specifying the minimum number of genes that a geneset must contain to be kept.
#' @param maxGenesetSize An integer specifying the maximum number of genes that a geneset can contain to be kept.
#'
#' @return A list (like the original GMT object) but only containing the genesets that met the size criteria.
#' @examples
#' # Consider a GMT object
#' gmt <- list(list(pathway = "geneset1", genes = c("gene1", "gene2", "gene3")),
#'              list(pathway = "geneset2", genes = c("gene4", "gene5")),
#'              list(pathway = "geneset3", genes = c("gene6", "gene7", "gene8", "gene9", "gene10")))
#'
#' # Filter the GMT object
#' filter_gmt_by_size(gmt, minGenesetSize = 3, maxGenesetSize = 4)
#' # This will return the first geneset only, as it is the only one with size between 3 and 4
#'
#' @export
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

#' Filter Gene Sets in a GMT File by Specific IDs
#'
#' This function filters gene sets contained within a GMT file based on a given vector of IDs.
#' It checks each gene set for a matching `external_id` and retains only those that match any of the provided IDs.
#'
#' @param gmt A list representing the GMT file's data.
#'           Each element of the list is expected to be a gene set, which is itself a list with an `external_id` field.
#' @param vector_id A vector of IDs (as character strings or numeric values).
#'                  Gene sets whose `external_id` matches any of these IDs will be retained in the output.
#'
#' @return A list of gene sets from the GMT file that match the specified IDs.
#'         This list will only contain elements where the `external_id` matched at least one of the IDs in `vector_id`.
#'         If no matches are found, an empty list is returned.
#'
#' @examples
#' # Assuming `gmt_data` is a list representing a GMT file
#' # and you want to filter for gene sets with external IDs 101 and 102:
#' filtered_data <- filter_gmt_by_id(gmt_data, c(101, 102))
#'
#' @note This function does not handle cases where `external_id` is missing in the gene sets.
#'       It is assumed that each gene set in the GMT has a valid `external_id` field.
#'
filter_gmt_by_id <- function(gmt, vector_id) {
  # Create an empty list to store filtered gene sets
  filtered_gmt <- list()

  # Loop over all gene sets in the GMT file
  for (i in seq_along(gmt)) {
    # Retrieve the external_id of the gene set
    external_id <- gmt[[i]]$external_id

    # Check if the external_id is in the provided vector_id
    if (external_id %in% vector_id) {
      filtered_gmt[[i]] <- gmt[[i]]
    }
  }

  # Remove NULL elements from the list
  filtered_gmt <- filtered_gmt[!sapply(filtered_gmt, is.null)]

  return(filtered_gmt)
}

#' Load CPDB GMT File
#'
#' This function loads the GMT file corresponding to the desired identifier from the ConsensusPathDB.
#'
#' @param identifier A character string specifying the desired identifier.
#' Options include 'ENSEMBL', 'Entrez', 'RefSeq', 'Symbol', 'Metabolites'.
#'
#' @return The loaded GMT file as a data frame.
#' @examples
#' loadCPDB("ENSEMBL")
#' @export
loadCPDB <- function(identifier) {

  # Check if the input is a character string
  if (!is.character(identifier)) {
    stop("The identifier must be a character string.")
  }

  # Check if the input is one of the allowed options
  allowed_identifiers <- c("ENSEMBL", "Entrez", "RefSeq", "Symbol", "Metabolites")
  if (!(identifier %in% allowed_identifiers)) {
    stop(paste("Invalid identifier. Options include:", paste(allowed_identifiers, collapse = ", ")))
  }

  # Load the corresponding GMT file
  if (identifier == "ENSEMBL") {
    data("CPDB_pathways_genes_ensembl")
    return(CPDB_pathways_genes_ensembl)
  } else if (identifier == "Entrez") {
    data("CPDB_pathways_genes_entrez")
    return(CPDB_pathways_genes_entrez)
  } else if (identifier == "RefSeq") {
    data("CPDB_pathways_genes_refseq")
    return(CPDB_pathways_genes_refseq)
  } else if (identifier == "Symbol") {
    data("CPDB_pathways_genes_symbol")
    return(CPDB_pathways_genes_symbol)
  } else if (identifier == "Metabolites") {
    data("CPDB_pathways_metabolites")
    return(CPDB_pathways_metabolites)
  }
}
