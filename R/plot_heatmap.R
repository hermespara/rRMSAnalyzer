#' Plot heatmap for a SummarizedExperiment object.
#' @description This easy function will let you display an heatmap for any given column (count or c-score). You can add an additionnal layer of information with metadata columns.
#' @param ribo A SummarizedExperiment object.
#' @param color_col Vector of the metadata columns’ name used for coloring samples.
#' @param sample_colors Optional custom colors for sample annotations. Supply a
#' named character vector when \code{color_col} has length 1, or a named list
#' of named character vectors when multiple metadata columns are used.
#' @param only_annotated Use only annotated sites (default = TRUE).
#' @param sites A character vector of specific annotated site names to use for
#' the heatmap (e.g., \code{c("28S_Am1310", "18S_Am99")}). Cannot be used
#' together with \code{only_annotated = TRUE}.
#' @param title Title to display on the plot. "default" for default title.
#' @param cutree_rows number of clusters the rows are divided into, based on the hierarchical clustering (using cutree).
#' @param cutree_cols number of clusters the columns are divided into, based on the hierarchical clustering (using cutree).
#' @param cluster_columns Column clustering specification passed to
#' \code{ComplexHeatmap::Heatmap()}. Can be a logical, an \code{hclust}
#' object, or a \code{dendrogram}. Default is \code{TRUE}.
#' @param sample_order Optional character vector specifying the exact order of
#' samples to display in the heatmap columns.
#' @param row_names_gp Font size for row names, passed to
#' \code{grid::gpar(fontsize = row_names_gp)}. Default is \code{6}.
#' @param clustering_distance Distance metric used for hierarchical clustering
#' of both rows and columns. Default is \code{"manhattan"}.
#' @param clustering_method Linkage method used for hierarchical clustering of
#' both rows and columns. Default is \code{"ward.D2"}.
#' @param ... Pheatmap’s parameters
#' @return A ggplot object of a heatmap. See ComplexHeatmap doc for more details
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # data("human_methylated")
#' # ribo_toy <- rename_rna(ribo_toy)
#' # ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' # plot_heatmap(ribo_toy,  color_col = c("run","condition"), only_annotated=TRUE)
plot_heatmap <- function(ribo, color_col = NULL, only_annotated = FALSE,
                         title = "default",
                         sites = NULL, cutree_rows = 4, cutree_cols = 2,
                         cluster_columns = TRUE,
                         sample_order = NULL,
                         row_names_gp = 6,
                         clustering_distance = "manhattan",
                         clustering_method = "ward.D2",
                         sample_colors = NULL, ...) {
  check_is_se(ribo)
  check_type(only_annotated, "logical", "only_annotated", length = 1)
  check_type(sites, "character", "sites")
  check_type(cutree_rows, "numeric", "cutree_rows", length = 1)
  check_type(cutree_cols, "numeric", "cutree_cols", length = 1)
  if (!is.logical(cluster_columns) &&
    !inherits(cluster_columns, "hclust") &&
    !inherits(cluster_columns, "dendrogram")) {
    cli::cli_abort(
      "{.arg cluster_columns} must be a logical, an {.cls hclust} object, or a {.cls dendrogram}."
    )
  }
  if (is.logical(cluster_columns) && length(cluster_columns) != 1) {
    cli::cli_abort("{.arg cluster_columns} must have length 1 when supplied as a logical.")
  }
  check_type(sample_order, "character", "sample_order")
  check_type(row_names_gp, "numeric", "row_names_gp", length = 1)
  validate_clustering_args(clustering_distance, clustering_method)
  check_sample(ribo, sample_order)

  check_metadata(ribo, color_col)
  matrix <- extract_data(ribo, "cscore",
    position_to_rownames = TRUE,
    only_annotated = only_annotated,
    sites = sites
  )

  .plot_heatmap(matrix, SummarizedExperiment::colData(ribo),
    color_col = color_col,
    most_variant = FALSE, title = title, cutree_rows = cutree_rows,
    cutree_cols = cutree_cols,
    cluster_columns = cluster_columns,
    sample_order = sample_order,
    row_names_gp = row_names_gp,
    clustering_distance = clustering_distance,
    clustering_method = clustering_method,
    sample_colors = sample_colors, ...
  )
}


#' Internal function to plot heatmap.
#'
#' @param cscore_matrix Sites x Samples C-score matrix (output of extract_data()).
#' @param metadata metadata for samples in cscore_matrix
#' @inheritParams plot_heatmap
#' @param most_variant select only the most variant positions (cannot be used from plot_heatmap())
#'
#' @return ComplexHeatmap heatmap
#' @keywords internal
.plot_heatmap <- function(cscore_matrix = NULL, metadata = NULL,
                          color_col = NULL, most_variant = FALSE,
                          title = "default", cutree_rows,
                          cutree_cols, cluster_columns = TRUE,
                          sample_order = NULL,
                          row_names_gp = 6,
                          clustering_distance = "manhattan",
                          clustering_method = "ward.D2",
                          sample_colors = NULL, ...) {
  heat_colors <- grDevices::hcl.colors(7, "inferno")
  plot_title <- resolve_plot_text(
    title,
    default_value = "Heatmap of C-score values",
    arg_name = "title"
  )


  column_ha <- build_heatmap_annotation(
    metadata = metadata,
    color_col = color_col,
    sample_colors = sample_colors
  )

  cscore_matrix <- stats::na.omit(cscore_matrix)
  cscore_matrix <- as.matrix(cscore_matrix)
  if (!is.null(sample_order)) {
    cscore_matrix <- cscore_matrix[, sample_order, drop = FALSE]
    metadata <- metadata[sample_order, , drop = FALSE]
    if (!is.null(column_ha)) {
      column_ha <- build_heatmap_annotation(
        metadata = metadata,
        color_col = color_col,
        sample_colors = sample_colors
      )
    }
  }
  ComplexHeatmap::Heatmap(cscore_matrix,
    col = heat_colors, name = "C-score",
    row_title = "Position", column_title = plot_title,
    column_title_side = "bottom",
    cluster_rows = TRUE, cluster_columns = cluster_columns,
    clustering_distance_columns = clustering_distance,
    clustering_distance_rows = clustering_distance,
    clustering_method_columns = clustering_method,
    clustering_method_rows = clustering_method,
    row_split = cutree_rows,
    column_split = if (cluster_columns) cutree_cols else NULL,
    top_annotation = column_ha,
    row_names_gp = grid::gpar(fontsize = row_names_gp),
    ...
  )
}
