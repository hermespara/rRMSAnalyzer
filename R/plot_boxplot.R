#' Plot a boxplot of a SummarizedExperiment counts.
#' @description This plot is useful to check if the samples are alike in their
#' raw counts.
#' @param ribo A SummarizedExperiment object.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param sample_colors Optional named character vector mapping metadata values
#' in \code{color_col} to colors.
#' @param outlier Show boxplot outlier values.
#' @param horizontal Show boxplot horizontally.
#' @return A ggplot object.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # boxplot_count(ribo_toy,"run")
#'
plot_boxplot_count <- function(ribo, color_col = NA,
                               outlier = TRUE, horizontal = FALSE,
                               sample_colors = NULL) {
  check_is_se(ribo)
  check_type(outlier, "logical", "outlier", length = 1)
  check_type(horizontal, "logical", "horizontal", length = 1)
  if (!is.na(color_col)) {
    validate_plot_colors(ribo, color_col, sample_colors)
  }

  ribo_matrix <- extract_data(ribo, "count", position_to_rownames = TRUE)

  return(.plot_boxplot_samples(ribo_matrix, "count",
    SummarizedExperiment::colData(ribo), color_col, outlier,
    sample_colors = sample_colors,
    horizontal = horizontal
  ))
}


#' Plot boxplot representing the C-score values of all samples for each
#' individual annotated site.
#' Sites are sorted by their median.
#' @inheritParams plot_boxplot_count
#' @param sort_by Sort sites by median ("median", default) by variance ("var")
#'  or IQR ("iqr").
#' @return a ggplot geom_boxplot
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy, human_methylated)
#' boxplot_cscores(ribo_toy)
#'
plot_boxplot_cscores <- function(ribo, outlier = TRUE, sort_by = c("median", "iqr", "var")[1], horizontal = FALSE) {
  check_is_se(ribo)
  check_type(outlier, "logical", "outlier", length = 1)
  check_in_set(tolower(sort_by), c("median", "iqr", "var"), "sort_by")
  check_type(horizontal, "logical", "horizontal", length = 1)

  ribo_m <- extract_data(ribo, only_annotated = TRUE)

  if (nrow(ribo_m) == 0) {
    cli::cli_abort("No annotated site found. Please use {.fn annotate_site} before calling this function.")
  }

  return(.plot_boxplot_sites(ribo_m,
    values_to_plot = "cscore", outlier = TRUE,
    sort_by = sort_by,
    horizontal = horizontal
  ))
}

#' Internal function for boxplot_cscore
#'
#' @param matrix Sites x Samples C-score/count matrix (output of extract_data()).
#' @param values_to_plot Value to display in plot.
#' @param outlier Show boxplot outlier values.
#' @param horizontal Show boxplot horizontally.
#'
#' @return A ggplot object.
#' @keywords internal
#'
.plot_boxplot_sites <- function(matrix, values_to_plot, outlier, sort_by, horizontal) {
  site <- cscore <- NULL
  id_vars <- "site"

  if (tolower(sort_by) == "iqr") {
    matrix_melted <- get_IQR(matrix)
  } else if (tolower(sort_by) == "median") {
    matrix_melted <- reshape2::melt(matrix,
      id.vars = id_vars,
      value.name = values_to_plot
    )
  } else if (tolower(sort_by) == "var") {
    matrix_melted <- get_IQR(matrix, "var")
  } else {
    stop(
      "Choose either \"median\", \"iqr\" or \"var\" for sort_by param.\n  \"",
      sort_by, "\" is not a valid argument."
    )
  }


  shape_outlier <- NA
  if (outlier) {
    shape_outlier <- 19
  }

  if ((tolower(sort_by) %in% c("iqr", "var"))) {
    method <- ifelse(tolower(sort_by) == "iqr", "IQR", "variance")
    p <- ggplot(matrix_melted, aes(x = site, y = cscore)) +
      geom_boxplot() +
      geom_boxplot() +
      theme_rRMSAnalyzer() +
      theme(
        axis.text.x = element_text(
          angle = 90, vjust = 0.5,
          hjust = 1
        ),
        legend.position = "top"
      ) +
      labs(
        x = paste0("Site (sorted by ", method, ")"),
        y = "C-score"
      )
  } else {
    matrix_melted <- reshape2::melt(matrix,
      id.vars = id_vars,
      value.name = values_to_plot
    )
    p <- ggplot2::ggplot(
      matrix_melted,
      ggplot2::aes(
        x = stats::reorder(site,
          !!rlang::sym(values_to_plot),
          na.rm = TRUE
        ),
        y = !!rlang::sym(values_to_plot)
      )
    ) +
      ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
      theme_rRMSAnalyzer() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5,
        hjust = 1
      )) +
      ggplot2::xlab("Site (sorted by median)") +
      ggplot2::ylab("C-score")
  }

  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  }

  return(p)
}

#' (internal) plot boxplot for a given matrix of values
#'
#' @param matrix Sites x Samples C-score/count matrix (output of extract_data()).
#' @param metadata Metadata of samples in matrix.
#' @param color_col Name of the column in the metadata used for coloring samples.
#' @param outlier Show boxplot outlier values.
#' @param values_col_name Name of the column containing the value (either count
#' or cscore).
#' @param horizontal Show boxplot horizontally.
#'
#' @return ggplot boxplot
#' @keywords internal
#'
.plot_boxplot_samples <- function(matrix, values_col_name, metadata,
                                  color_col = NA, outlier, horizontal,
                                  sample_colors = NULL) {
  Sample <- NULL
  id_vars <- "Sample"
  matrix <- log10(matrix)
  matrix_inv <- as.data.frame(t(matrix))
  # matrix_inv <- tibble::rownames_to_column(matrix_inv, "Sample")
  matrix_inv["Sample"] <- rownames(matrix_inv)

  if (!is.na(color_col)) {
    matrix_inv <- cbind(matrix_inv, metadata[color_col])
    id_vars <- c(color_col, "Sample")
  }

  matrix_melted <- reshape2::melt(matrix_inv,
    id.vars = id_vars,
    value.name = values_col_name
  )

  matrix_melted[["Sample"]] <- factor(matrix_melted[["Sample"]],
    levels = unique(matrix_melted[["Sample"]])
  )
  shape_outlier <- NA

  if (outlier) {
    shape_outlier <- 19
  }

  if (is.na(color_col)) {
    p <- ggplot2::ggplot(
      matrix_melted,
      ggplot2::aes(x = Sample, y = !!sym(values_col_name))
    )
  } else {
    p <- ggplot2::ggplot(
      matrix_melted,
      ggplot2::aes(
        x = Sample, y = !!sym(values_col_name),
        fill = !!sym(color_col)
      )
    )
  }
  p <- p + ggplot2::geom_boxplot(outlier.shape = shape_outlier) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, vjust = 0.5,
      hjust = 1
    ))

  if (!is.na(color_col)) {
    p <- p + plot_group_scales(sample_colors = sample_colors, color = FALSE, fill = TRUE)
  }


  p <- p + ggplot2::geom_hline(yintercept = 2, colour = "#56B4E9", linetype = "dashed")
  if (horizontal) {
    p <- p + ggplot2::coord_flip()
  }

  return(p)
}
