# ==================================================

#' Plot a correlation heatmap from a SummarizedExperiment object.
#'
#' Shows the correlation **distance** between samples.
#' @md
#' @inheritParams plot_heatmap
#' @param values_col Name of the column containing the value (either count or cscore).
#' @param title Title to display on the heatmap. Use \code{"default"} for the
#' standard title.
#' @return ComplexHeatmap object
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # plot_heatmap_corr(ribo_toy,"count","condition")
#'
plot_heatmap_corr <- function(ribo, values_col, color_col = NULL,
                              title = "default",
                              clustering_distance = "pearson",
                              clustering_method = "ward.D2",
                              sample_colors = NULL) {
  check_is_se(ribo)
  check_type(values_col, "character", "values_col", length = 1)
  validate_clustering_args(clustering_distance, clustering_method)
  matrix <- extract_data(ribo, values_col, position_to_rownames = TRUE)
  if (!is.null(color_col)) {
    check_metadata(ribo, color_col)
  }
  .plot_heatmap_corr(
    matrix,
    SummarizedExperiment::colData(ribo),
    color_col = color_col,
    title = title,
    clustering_distance = clustering_distance,
    clustering_method = clustering_method,
    sample_colors = sample_colors
  )
}
#' Internal function of plot_heatmap_corr.
#'
#' @param cscore_matrix  Sites x Samples C-score matrix (output of extract_data()).
#' @param metadata Metadata of samples in matrix
#' @param color_col Vector of the metadata columns’ name used for coloring samples.
#'
#' @return ComplexHeatmap heatmap
#' @keywords internal
#'
.plot_heatmap_corr <- function(cscore_matrix, metadata,
                               color_col, title = "default",
                               clustering_distance = "pearson",
                               clustering_method = "ward.D2",
                               sample_colors = NULL) {
  column_ha <- build_heatmap_annotation(
    metadata = metadata,
    color_col = color_col,
    sample_colors = sample_colors,
    na_col = "red"
  )
  plot_title <- resolve_plot_text(
    title,
    default_value = "Pearson correlation between samples",
    arg_name = "title"
  )


  corr_matrix <- stats::cor(cscore_matrix, use = "complete.obs")
  pearson_color <- colorRamp2::colorRamp2(c(0, 0.5, 1), c("red", "white", "blue"))

  ComplexHeatmap::Heatmap(corr_matrix,
    col = pearson_color, name = "Pearson correlation",
    row_title = "Sample", column_title = plot_title,
    column_title_side = "top",
    cluster_rows = FALSE, cluster_columns = TRUE,
    clustering_distance_columns = clustering_distance,
    heatmap_legend_param = list(
      legend_direction = "horizontal"
    ),
    clustering_method_columns = clustering_method,
    column_split = 3,
    top_annotation = column_ha,
    row_names_gp = grid::gpar(fontsize = 6)
  )
}
