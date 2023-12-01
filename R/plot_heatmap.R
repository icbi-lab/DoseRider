#' Create Complex Heatmap for Gene Sets
#'
#' This function creates a heatmap using the `ComplexHeatmap` package, where each row represents a gene set
#' and each column represents the average expression of genes within that gene set for each dose.
#'
#' @param dose_rider_results A list containing the results of the DoseRider analysis for each gene set.
#' Each element of the list is a sublist with various metrics and the raw values (expression data) for a gene set.
#' @param dose_col A character string specifying the name of the dose column in the raw expression data.
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
dose_response_heatmap <- function(dose_rider_results, dose_col, top = 15) {
  # Initialize an empty matrix to store the average expressions
  heatmap_data <- list()

  #Top pathways in function of P-Value
  dose_rider_df <- as.data.frame(dose_rider_results)
  dose_rider_df <- dose_rider_df[order(dose_rider_df$Adjusted_Cubic_P_Value),]

  # Gene set names
  gene_set_names <- dose_rider_df$Geneset[1:top]

  for (gene_set_name in gene_set_names) {
    res_geneset <- dose_rider_results[[gene_set_name]]
    gene_set_name <- str_replace(gene_set_name, " - Homo sapiens \\(human\\)", "")
    # Extract the raw values
    raw_values <- res_geneset$Raw_Values[[1]]

    # Calculate the average expression for each dose
    avg_expression <- aggregate(counts ~ get(dose_col, raw_values), data = raw_values, FUN = mean)
    rownames(avg_expression) <- avg_expression$`get(dose_col, raw_values)`
    # Store the average expressions with gene set name as column name
    heatmap_data[[gene_set_name]] <- setNames(avg_expression[["counts"]], avg_expression[[dose_col]])
  }

  # Combine the list into a matrix
  heatmap_matrix <- do.call(cbind, heatmap_data)

  # Calculate Z-scores across all doses for each gene set
  z_score_matrix <- t(apply(t(heatmap_matrix), 1, function(x) (x - mean(x)) / sd(x)))
  rownames(z_score_matrix) <- unlist(lapply(colnames(heatmap_matrix),function(x){str_wrap(x,width = 35)}))
  colnames(z_score_matrix) <- rownames(avg_expression)

  # Create the heatmap
  col_fun <- colorRamp2(c(-2, -1, 0, 1, 2), c("blue", "lightblue", "white", "lightcoral", "red"))

  ha <- Heatmap(z_score_matrix,
                name = "Z-Score",
                column_title = "Dose",
                column_title_side = "bottom",
                row_title = "Gene Set",
                column_gap = unit(2, "mm"),
                border_gp = grid::gpar(col = "black", lty = 1),
                rect_gp = grid::gpar(col = "black", lwd = 1),
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_gp = gpar(fontsize = 10, fontface = "bold", just = "center"),
                column_names_rot = 0,
                cluster_columns = FALSE,
                show_row_dend = FALSE,
                col = col_fun)


  return(ha)
  }

