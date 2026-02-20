#' Principal component analysis of a SummarizedExperiment object
#'
#'
#' @param ribo A SummarizedExperiment object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param axes Two-element vector indicating which pair of principal components
#' to show.
#' @param only_annotated If TRUE, use only annotated sites to plot PCA.
#' @param sites A character vector of specific annotated site names to use for
#' PCA (e.g., \code{c("28S_Am1310", "18S_Am99")}). Cannot be used together with
#' \code{only_annotated = TRUE}.
#' @param title Title to display on the plot. 'default' for default title.
#' @param subtitle Subtitle to display on the plot. 'samples' for number of
#' samples. 'none' for no subtitle.
#' @param draw_ellipses If TRUE, draw ellipses around groups.
#' @param object_only Return directly the full dudi.pca object, without
#' generating the plot.
#' @return A ggplot or a dudi.pca object if object_only is set to True.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # plot_pca(ribo_toy,'run')
#' plot_pca(ribo_toy, "run", draw_ellipses = TRUE)
plot_pca <- function(ribo, color_col = NULL, axes = c(1, 2),
                     only_annotated = FALSE, sites = NULL,
                     title = "default",
                     subtitle = "samples", draw_ellipses = FALSE,
                     object_only = FALSE) {
  if (missing(ribo)) {
    cli::cli_abort(c(
      "Missing argument : {.var ribo}",
      "i" = "{.var ribo} only accepts a SummarizedExperiment object."
    ))
  }
  if (!inherits(ribo, "SummarizedExperiment")) {
    cli::cli_abort(c(
      "{.var ribo} must be a SummarizedExperiment object",
      "x" = "You've supplied a {.cls {class(ribo)}}."
    ))
  }
  if (isFALSE(ribo@metadata$has_cscore)) {
    cli::cli_abort(c(
      "No C-score found in the object supplied in {.var ribo}!",
      "i" = "You can compute C-scores using compute_cscore function."
    ))
  }

  if (!is.null(color_col)) {
    check_metadata(ribo, color_col)
  }

  check_type(sites, "character", "sites")

  pca_matrix <- extract_data(ribo, "cscore",
    position_to_rownames = TRUE,
    only_annotated = only_annotated,
    sites = sites
  )

  pca_calculated <- .calculate_pca(pca_matrix)

  if (object_only) {
    return(pca_calculated)
  }

  facto_pca <- .plot_pca(pca_calculated, SummarizedExperiment::colData(ribo), color_col,
    axes = axes, title = title, subtitle = subtitle, draw_ellipses = draw_ellipses
  )
  return(facto_pca)
}

#' Compute PCA from a c-score matrix
#'
#' @keywords internal
#'
#' @param cscore.matrix matrix of c-score extracted from a SummarizedExperiment with ' \code{\link{extract_data}}
#' @return dudi.pca object
#'
.calculate_pca <- function(cscore.matrix = NULL) {
  pca.res <- ade4::dudi.pca(t(cscore.matrix[stats::complete.cases(cscore.matrix), ]),
    scannf = FALSE,
    nf = 5
  )

  return(pca.res)
}

#' Plot a dudi.pca object using ggplot2
#'
#' @keywords internal
#'
#' @param dudi.pca a dudi.pca object generated with Ade4
#' @param metadata metadata table from SummarizedExperiment
#' @param draw_ellipses If TRUE, draw ellipses around groups.
#' @inheritParams plot_pca
#'
#' @return a ggplot
#'
.plot_pca <- function(dudi.pca = NULL, metadata = NULL,
                      color_col = NULL, axes = axes, title = "default",
                      subtitle = "samples", draw_ellipses = FALSE) {
  # Prepare data for plotting
  df_pca <- data.frame(dudi.pca$li)
  colnames(df_pca) <- paste0("Axis", seq_len(ncol(df_pca)))

  if (!is.null(metadata)) {
    # Ensure metadata matches PCA rows if needed, assuming order is preserved from input
    # dudi.pca$li rows correspond to input matrix columns (samples)
    df_pca <- cbind(df_pca, metadata)
  }

  if (is.null(color_col)) {
    color_column <- "none"
  } else {
    color_column <- metadata[, color_col]
  }

  # Title
  if (title == "default") {
    plot_title <- "Principal Component Analysis from C-score data"
  } else {
    plot_title <- title
  }

  # Subtitle
  if (subtitle == "samples") {
    plot_subtitle <- paste(
      nrow(df_pca),
      "samples and", ncol(dudi.pca$tab),
      "positions"
    )
  } else if (subtitle == "none") {
    plot_subtitle <- ggplot2::waiver()
  } else {
    plot_subtitle <- subtitle
  }

  # Calculate variance explained
  eig_percent <- round(dudi.pca$eig / sum(dudi.pca$eig) * 100, 1)

  # Define axes columns
  x_axis <- paste0("Axis", axes[1])
  y_axis <- paste0("Axis", axes[2])

  # Base plot
  p <- ggplot2::ggplot(df_pca, ggplot2::aes(x = .data[[x_axis]], y = .data[[y_axis]]))

  if (!identical(color_column, "none")) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[color_col]]), size = 3) +
      scale_color_rRMSAnalyzer() +
      scale_fill_rRMSAnalyzer()

    if (draw_ellipses) {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(color = .data[[color_col]], fill = .data[[color_col]]),
        geom = "polygon", alpha = 0.1, show.legend = FALSE
      )
    }
  } else {
    p <- p + ggplot2::geom_point(size = 3)
  }

  # Add text labels (simple repel effect or just text)
  # Since we removed factoextra/ggrepel to save deps, we use geom_text with check_overlap or vjust
  # Or we can check if ggrepel is available. Since we want to be lightweight, we use geom_text.
  p <- p + ggplot2::geom_text(ggplot2::aes(label = rownames(df_pca)), vjust = -1, size = 4, check_overlap = FALSE)

  p <- p +
    theme_rRMSAnalyzer() +
    ggplot2::labs(title = plot_title, subtitle = plot_subtitle) +
    ggplot2::xlab(paste0("PC", axes[1], ": ", eig_percent[axes[1]], "%")) +
    ggplot2::ylab(paste0("PC", axes[2], ": ", eig_percent[axes[2]], "%"))

  return(p)
}
