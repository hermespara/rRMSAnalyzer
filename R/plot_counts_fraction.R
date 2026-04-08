#' Plot count distribution among RNAs for each sample
#'
#' @param ribo A SummarizedExperiment object.
#' @param title Title to display on the plot. Use \code{"default"} for the
#' standard title.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # plot_counts_fraction(ribo_toy)
plot_counts_fraction <- function(ribo, title = "default") {
  check_is_se(ribo)
  plot_title <- resolve_plot_text(
    title,
    default_value = "Counts distribution per sample",
    arg_name = "title"
  )
  # NSE fix
  samplesid <- rna <- count <- sum.counts <- NULL

  # Get counts matrix
  mat <- SummarizedExperiment::assay(ribo, "counts")

  # Get RNA info
  rd <- SummarizedExperiment::rowData(ribo) # assumed to align with matrix rows

  # Pivot to long format for aggregation
  # Need to include RNA info
  df_wide <- as.data.frame(mat)
  df_wide$rna <- rd$rna

  df_long <- tidyr::pivot_longer(df_wide, cols = -rna, names_to = "samplesid", values_to = "count")

  # Aggregation
  all_data_sums <- df_long |>
    dplyr::group_by(samplesid, rna) |>
    dplyr::summarise(
      sum.counts = sum(count, na.rm = TRUE),
      # n0 logic from original: length(which(count < 5))
      n0 = length(which(count < 5)),
      .groups = "drop"
    )

  rna_names_df <- S4Vectors::metadata(ribo)$rna_names

  if (!is.null(rna_names_df)) {
    all_data_sums$rna <- factor(all_data_sums$rna, levels = rna_names_df$current_name)
  }

  all_data_sums$samplesid <- factor(all_data_sums$samplesid, levels = colnames(ribo))

  counts_plot <- ggplot2::ggplot(all_data_sums, aes(x = samplesid, y = sum.counts, fill = rna))

  counts_plot <- counts_plot + ggplot2::geom_bar(
    stat = "identity",
    position = "fill"
  ) +
    scale_fill_rRMSAnalyzer() +
    theme_rRMSAnalyzer() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 1)) +
    labs(
      title = plot_title,
      x = "Sample",
      y = "Counts fraction"
    )
  return(counts_plot)
}
