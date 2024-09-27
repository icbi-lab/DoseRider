#' Create Complex Heatmap for Gene Sets
#'
#' This function creates a heatmap using the `ComplexHeatmap` package, where each row represents a gene set
#' and each column represents the average expression of genes within that gene set for each dose.
#'
#' @param dose_rider_results A list containing the results of the DoseRider analysis for each gene set.
#' Each element of the list is a sublist with various metrics and the raw values (expression data) for a gene set.
#' @param dose_col A character string specifying the name of the dose column in the raw expression data.
#' @param dose_unit A character string specifying the dose units to plot in the column title..
#' @param top An integer specifying the number of top gene sets to include in the heatmap. Default is 15.
#' @param order_column A character string specifying the column to use for ordering gene sets in the heatmap.
#' @param fontsize Integer for fontsize. Default 6
#'
#' @return An object of class `Heatmap` representing the constructed heatmap.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom stringr str_replace str_wrap
#'
#' @examples
#' \dontrun{
#'   # Assuming `dose_rider_results` is your DoseRider analysis result
#'   heatmap_plot <- create_complex_heatmap(dose_rider_results, "dose")
#'   draw(heatmap_plot) # To display the heatmap
#' }
#'
#' @export
dose_response_heatmap <- function(dose_rider_results, dose_col = "Dose", dose_unit = "μM", top = 15, order_column = "NegLogPValue", decreasing = FALSE, fontsize = 6) {
  # Initialize an empty matrix to store the average expressions
  heatmap_data <- list()

  #Top pathways in function of P-Value
  # Select top gene sets
  top_gene_sets <- get_top_genesets(dose_rider_results = dose_rider_results,
                                    top = top,
                                    decreasing = decreasing,
                                    order_column = order_column )

  for (gene_set_name in top_gene_sets) {
    res_geneset <- dose_rider_results[[gene_set_name]]
    gene_set_name <- str_replace(gene_set_name, " - Homo sapiens \\(human\\)", "")
    # Extract the raw values
    raw_values <- res_geneset$Raw_Values[[1]]

    # Calculate the average expression for each dose
    if ("predictions" %in% colnames(raw_values)) {
      avg_expression <- aggregate(predictions ~ get(dose_col, raw_values), data = raw_values, FUN = mean)
      rownames(avg_expression) <- avg_expression$`get(dose_col, raw_values)`
      gene_set_name <- gsub("_", " ", gene_set_name)
      # Store the average expressions with gene set name as column name
      heatmap_data[[gene_set_name]] <- setNames(avg_expression[["predictions"]], avg_expression[[dose_col]])
    }
  }

  # Combine the list into a matrix
  heatmap_matrix <- do.call(cbind, heatmap_data)

  # Calculate Z-scores across all doses for each gene set
  z_score_matrix <- t(apply(t(heatmap_matrix), 1, function(x) (x - mean(x)) / sd(x)))
  rownames(z_score_matrix) <- unlist(lapply(colnames(heatmap_matrix), function(x) {str_wrap(x, width = 20)}))
  colnames(z_score_matrix) <- round(as.numeric(rownames(avg_expression)), 3)

  # Create the heatmap
  col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("blue", "lightblue", "white", "lightcoral", "red"))

  ha <- Heatmap(z_score_matrix,
                name = "Z-Score",
                column_title = paste0(dose_col, " (",dose_unit,")"),
                column_title_side = "bottom",
                column_gap = unit(2, "mm"),
                border_gp = grid::gpar(col = "black", lty = 1),
                rect_gp = grid::gpar(col = "black", lwd = 1),
                row_names_gp = gpar(fontsize = fontsize),  # Matches axis.text size in theme
                column_names_gp = gpar(fontsize = fontsize, just = "center"),  # Matches axis.text size in theme
                column_title_gp = gpar(fontsize = fontsize + 2, just = "center"),
                column_names_rot = 45,
                cluster_columns = FALSE,
                show_row_dend = FALSE,
                col = col_fun,
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = fontsize -2),  # Matches legend.title size in theme
                  labels_gp = gpar(fontsize = fontsize-2),  # Matches legend.text size in theme
                  legend_direction = "horizontal",
                  legend_position = "bottom"
                ))

  return(ha)
}

#' Create Heatmap for Individual Genes within a Gene Set
#'
#' This function creates a heatmap for a specified gene set from DoseRider results. Each row in the heatmap
#' represents an individual gene, and each column represents different doses, showing gene expression values.
#'
#' @param dose_rider_results A list containing the results of the DoseRider analysis for each gene set.
#' @param gene_set_name A character string specifying the name of the gene set to be visualized.
#' @param dose_col A character string specifying the name of the dose column in the raw expression data.
#' @param dose_unit A character string specifying the dose units to plot in the column title..
#' @param fontsize Integer for fontsize. Default 6
#'
#' @return An object of class `Heatmap` representing the constructed heatmap for the specified gene set.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#' @importFrom reshape2 dcast
#'
#' @examples
#' \dontrun{
#'   # Assuming `dose_rider_results` is your DoseRider analysis result
#'   heatmap_plot <- create_gene_heatmap(dose_rider_results, "Gene Set Name", "dose")
#'   draw(heatmap_plot) # To display the heatmap
#' }
#'
#' @export
create_gene_heatmap <- function(dose_rider_results, gene_set_name, dose_col,  dose_unit = "μM", fontsize = 6) {
  if (!gene_set_name %in% names(dose_rider_results)) {
    stop("Specified gene set name not found in the results.")
  }

  # Extract the raw values for the specified gene set
  res_geneset <- dose_rider_results[[gene_set_name]]
  raw_values <- res_geneset$Raw_Values[[1]]

  # Calculate the average expression for each dose
  if ("predictions" %in% colnames(raw_values)){
    # Calculate the average expression for each gene at each dose
    avg_expression <- aggregate(predictions ~ get(dose_col, raw_values) + gene,
                                data = raw_values,
                                FUN = mean)
    colnames(avg_expression) <- c(dose_col, "gene", "predictions")

    # Pivot the data to have genes as rows and doses as columns
    heatmap_data <- dcast(avg_expression, gene ~ get(dose_col, avg_expression),
                          value.var = "predictions")
    rownames(heatmap_data) <- heatmap_data$gene
    heatmap_data$gene <- NULL
  }



  # Calculate Z-scores for standardization
  z_score_matrix <- t(apply(heatmap_data, 1, scale))
  colnames(z_score_matrix) <- colnames(heatmap_data)
  # Define color function for the heatmap

  col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("blue", "lightblue", "white", "lightcoral", "red"))
  names(custom_palette) <- as.character(1:length(custom_palette))

  if("ClusterAssignments" %in% names(res_geneset)){
    row_ha = HeatmapAnnotation(Cluster = res_geneset$ClusterAssignments,
                               col = list(Cluster = custom_palette), which = "row",
                               show_annotation_name = FALSE,
                               show_legend = F)
  } else {
    row_ha = rowAnnotation(foo = anno_empty(border = FALSE,
                                            width = unit(0.1, "mm")))
  }

  # Create the heatmap
  ha <- Heatmap(z_score_matrix,
                name = "Z-Score",
                column_title = paste0(dose_col, " (",dose_unit,")"),
                column_title_side = "bottom",
                column_gap = unit(2, "mm"),
                border_gp = grid::gpar(col = "black", lty = 1),
                rect_gp = grid::gpar(col = "black", lwd = 1),
                row_names_gp = gpar(fontsize = fontsize),  # Matches axis.text size in theme
                column_names_gp = gpar(fontsize = fontsize, just = "center"),  # Matches axis.text size in theme
                column_title_gp = gpar(fontsize = fontsize + 2, just = "center"),
                column_names_rot = 45,
                right_annotation = row_ha,
                cluster_columns = FALSE,
                show_row_dend = FALSE,
                col = col_fun,
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = fontsize),  # Matches legend.title size in theme
                  labels_gp = gpar(fontsize = fontsize),  # Matches legend.text size in theme
                  legend_direction = "horizontal",
                  legend_position = "bottom"
                ))

  return(ha)
}
