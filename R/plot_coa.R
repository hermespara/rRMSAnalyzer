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
                     subtitle = "default", object_only = FALSE) {
  check_is_se(ribo)
  check_type(only_annotated, "logical", "only_annotated", length = 1)
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
    axes = axes, title, subtitle
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
                      subtitle = "default") {
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
    # color_col is a vector of values here, but for aes mapping it's better to use column name if possible
    # but we extracted it as a vector.
    # To use scale_color_rRMSAnalyzer properly, we should map to the column in df_coa if it exists
    # We added metadata to df_coa, so we can use .data[[color_colname]]
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[color_colname]]), size = 2) +
      scale_color_rRMSAnalyzer() +
      ggplot2::labs(color = color_colname)
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
