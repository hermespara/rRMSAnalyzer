#' Correspondence analysis of a SummarizedExperiment' counts.
#'
#' @param ribo A SummarizedExperiment object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param axes Two-element vector indicating which pair of COA components
#' to show.
#' @param only_annotated If TRUE, use only annotated sites to plot COA.
#' @param title Title to display on the plot. 'default' for default title.
#' @param subtitle Subtitle to display on the plot. 'samples' for number of
#' samples. 'none' for no subtitle.
#' @param draw_ellipses If TRUE, draw ellipses around groups.
#' @param draw_centroids If TRUE, draw group centroids on the plot.
#' @param object_only Return directly the full dudi.coa object, without
#' generating the plot.
#' @return A ggplot or a dudi.coa object if object_only is set to True.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # plot_coa(ribo = ribo_toy, color_col = 'condition')
plot_coa <- function(ribo, color_col = NULL, axes = c(1, 2),
                     only_annotated = FALSE, title = "default",
                     subtitle = "default", draw_ellipses = FALSE,
                     draw_centroids = FALSE, object_only = FALSE) {
  check_is_se(ribo)
  check_type(only_annotated, "logical", "only_annotated", length = 1)
  check_type(draw_ellipses, "logical", "draw_ellipses", length = 1)
  check_type(draw_centroids, "logical", "draw_centroids", length = 1)
  check_type(object_only, "logical", "object_only", length = 1)
  if (!is.null(color_col)) {
    check_metadata(ribo, color_col)
  }

  coa_matrix <- extract_data(ribo, "counts",
    position_to_rownames = TRUE,
    only_annotated = only_annotated
  )

  coa_calculated <- .compute_coa(coa_matrix)

  if (object_only) {
    return(coa_calculated)
  }

  return(.plot_coa(coa_calculated, SummarizedExperiment::colData(ribo),
    color_col,
    axes = axes, title, subtitle,
    draw_ellipses = draw_ellipses, draw_centroids = draw_centroids
  ))
}

#' (internal) Correspondence Analysis computation function
#'
#' @param raw_counts A matrix of counts, as exported by extract_data()
#' @return A coa object
#' @keywords internal
#'
.compute_coa <- function(raw_counts = NULL) {
  res_coa <- ade4::dudi.coa(raw_counts[stats::complete.cases(raw_counts), ], scannf = FALSE, nf = 5)

  return(res_coa)
}


#' (internal) Plot a coa object
#'
#' @param dudi.coa A coa object
#' @param metadata SummarizedExperiment's metadata dataframe
#' @inheritParams plot_coa
#' @return A ggplot object containing the COA
#' @keywords internal
#'
#' (internal) Plot a coa object
#'
#' @param dudi.coa A coa object
#' @param metadata SummarizedExperiment's metadata dataframe
#' @inheritParams plot_coa
#' @return A ggplot object containing the COA
#' @keywords internal
#'
.plot_coa <- function(dudi.coa = NULL,
                      metadata = NULL, color_col = NULL,
                      axes = c(1, 2), title = "default",
                      subtitle = "default", draw_ellipses = FALSE,
                      draw_centroids = FALSE) {
  # Prepare data for plotting (samples are columns in COA of sites x samples)
  # dudi.coa on (sites x samples) -> $co are column coordinates (samples)
  df_coa <- data.frame(dudi.coa$co)
  colnames(df_coa) <- paste0("Axis", seq_len(ncol(df_coa)))

  # Ensure metadata matches COA columns (samples)
  # We assume the order is preserved.
  if (!is.null(metadata)) {
    df_coa <- cbind(df_coa, metadata)
  }

  if (is.null(color_col)) {
    color_colname <- color_col <- "Black"
  } else {
    color_colname <- color_col
    color_col <- metadata[, color_col]
  }

  if (title == "default") {
    title <- "Correspondance analysis from count data"
  }

  if (subtitle == "default") {
    subtitle <- paste(
      ncol(dudi.coa$tab), "samples and", nrow(dudi.coa$tab),
      "positions (all)"
    )
  }

  # Calculate variance explained (eigenvalues)
  eig_percent <- round(dudi.coa$eig / sum(dudi.coa$eig) * 100, 1)

  # Define axes columns
  x_axis <- paste0("Axis", axes[1])
  y_axis <- paste0("Axis", axes[2])

  # Base Plot
  p <- ggplot2::ggplot(df_coa, ggplot2::aes(x = .data[[x_axis]], y = .data[[y_axis]]))

  # Color handling
  if (!identical(color_colname, "Black")) {
    group_layers <- .group_ellipses(df_coa, x_axis, y_axis, color_colname)

    if (draw_ellipses && !is.null(group_layers$ellipses)) {
      p <- p + ggplot2::geom_polygon(
        data = group_layers$ellipses,
        ggplot2::aes(x = .data[["x"]], y = .data[["y"]], group = .data[["group"]],
          color = .data[["group"]], fill = .data[["group"]]
        ),
        inherit.aes = FALSE, alpha = 0.1, show.legend = FALSE
      )
    }

    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[color_colname]]), size = 2) +
      scale_color_rRMSAnalyzer() +
      scale_fill_rRMSAnalyzer() +
      ggplot2::labs(color = color_colname)

    if (draw_centroids && !is.null(group_layers$centroids)) {
      p <- p + ggplot2::geom_point(
        data = group_layers$centroids,
        ggplot2::aes(x = .data[["x"]], y = .data[["y"]], color = .data[["group"]]),
        inherit.aes = FALSE, shape = 4, stroke = 1.2, size = 4, show.legend = FALSE
      )
    }
  } else {
    p <- p + ggplot2::geom_point(color = "black", size = 2)
  }

  # Add text labels
  p <- p + ggplot2::geom_text(ggplot2::aes(label = rownames(df_coa)), vjust = -1, size = 4, check_overlap = FALSE)

  p <- p + theme_rRMSAnalyzer() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::xlab(paste0("Axis ", axes[1], ": ", eig_percent[axes[1]], "%")) +
    ggplot2::ylab(paste0("Axis ", axes[2], ": ", eig_percent[axes[2]], "%"))

  return(p)
}
