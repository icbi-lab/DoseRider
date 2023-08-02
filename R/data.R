#' #' CPDB Pathways and Genes (ENSEMBL)
#' #'
#' #' A dataset containing the ConsensusPathDB (CPDB) pathways and their associated genes
#' #' identified by ENSEMBL gene IDs.
#' #'
#' #' @format A data frame with the following columns:
#' #' \describe{
#' #'   \item{Pathway}{A character vector with the name of the pathway.}
#' #'   \item{Gene}{A character vector with the ENSEMBL gene IDs associated with the pathway.}
#' #' }
#' #' @source \url{http://consensuspathdb.org/}
#' "CPDB_pathways_genes_ensembl"

#' CPDB Pathways and Genes (Entrez)
#'
#' A dataset containing the ConsensusPathDB (CPDB) pathways and their associated genes
#' identified by Entrez gene IDs.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Pathway}{A character vector with the name of the pathway.}
#'   \item{Gene}{A character vector with the Entrez gene IDs associated with the pathway.}
#' }
#' @source \url{http://consensuspathdb.org/}
"CPDB_pathways_genes_entrez"

#' CPDB Pathways and Genes (RefSeq)
#'
#' A dataset containing the ConsensusPathDB (CPDB) pathways and their associated genes
#' identified by RefSeq gene IDs.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Pathway}{A character vector with the name of the pathway.}
#'   \item{Gene}{A character vector with the RefSeq gene IDs associated with the pathway.}
#' }
#' @source \url{http://consensuspathdb.org/}
"CPDB_pathways_genes_refseq"


#' CPDB Pathways and Metabolites
#'
#' A dataset containing the ConsensusPathDB (CPDB) pathways and their associated metabolites.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Pathway}{A character vector with the name of the pathway.}
#'   \item{Metabolite}{A character vector with the identifiers of metabolites associated with the pathway.}
#' }
#' @source \url{http://consensuspathdb.org/}
"CPDB_pathways_metabolites"

#' CPDB Pathways and Genes (Symbol)
#'
#' A dataset containing the ConsensusPathDB (CPDB) pathways and their associated genes
#' identified by Symbol gene IDs.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{Pathway}{A character vector with the name of the pathway.}
#'   \item{Gene}{A character vector with the Symbol gene IDs associated with the pathway.}
#' }
#' @source \url{http://consensuspathdb.org/}
"CPDB_pathways_genes_symbol"

#' Transcriptomic Data from Study PRJNA869442: BPA Alternatives in Cultured Breast Cancer Cells
#'
#' Overview:
#' Study PRJNA869442 is a transcriptomics study that aims to evaluate potential hazards and compare potencies of Bisphenol A (BPA) and 15 BPA alternative chemicals in cultured breast cancer cells (MCF-7). The study uses high-throughput transcriptomics to examine general toxicological effects and estrogen receptor alpha (ERα)-associated transcriptional changes in response to chemical exposures.
#'
#' Study Design:
#' - Cell Line: MCF-7 breast cancer cells
#' - Exposure: Cells were exposed to BPA and 15 BPA alternative chemicals at concentrations ranging from 0.0005 to 100 µM.
#' - Duration: Exposure period was 48 hours.
#' - Controls: The study includes technical controls (reference RNA, media only, cells with no treatment) and solvent controls (0.1% DMSO).
#' - Reference Chemicals: The study also includes two reference chemicals, estradiol, and dexamethasone.
#'
#' Data Format:
#' - Data Type: Gene expression data (RNA-Seq)
#' - Data Format: The data is provided as a SummarizedExperiment object or a matrix/data frame with metadata.
#' - Gene Sets: The data includes information about gene sets and their corresponding expression levels in response to different doses of BPA and its alternative chemicals.
#'
#' Preprocessing:
#' Before conducting the doseRider analysis, the transcriptomic data from PRJNA869442 has been preprocessed locally to ensure the data's quality and readiness for analysis.
#'
#' Objective:
#' The objective of this analysis is to utilize doseRider, a function in the doseRider package, to investigate potential non-linear dose-response relationships in the transcriptomic data. By using Generalized Additive Mixed Models (GAMMs), doseRider allows us to understand the relationship between dose-response and gene expression, providing insights into the biological effects and potencies of BPA and its alternative chemicals in cultured breast cancer cells.
#'
#' Data Source:
#' The original data for study PRJNA869442 can be accessed from the NCBI Sequence Read Archive (SRA). However, for this analysis, we have already preprocessed the data locally and will be using the preprocessed dataset.
#'
#' @format A SummarizedExperiment object or a matrix/data frame with metadata.
#' @examples
#' # Load the preprocessed dataset
#' data("PRJNA869442")
#'
#' # Check the structure of the dataset
#' str(PRJNA869442)
#'
#' # Perform analysis using doseRider
#' result <- doseRider(PRJNA869442)
#'
#' @seealso
#' \code{\link{doseRider}}
#'
#' @source The original data for study PRJNA869442 can be accessed from the NCBI Sequence Read Archive (SRA).
#'
"PRJNA869442"

