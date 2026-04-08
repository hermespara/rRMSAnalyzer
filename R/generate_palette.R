# Generate the col argument for ComplexHeatmap::HeatmapAnnotation
generate_palette <- function(metadata, cols_to_use, custom_colors = NULL) {
  annot_list <- list()
  default_palette <- .rRMSAnalyzer_palette()

  if (!is.null(custom_colors) && !is.list(custom_colors)) {
    if (length(cols_to_use) != 1) {
      cli::cli_abort(
        "{.arg sample_colors} must be a named list when {.arg color_col} contains multiple metadata columns."
      )
    }
    custom_colors <- stats::setNames(list(custom_colors), cols_to_use)
  } else if (is.list(custom_colors) && length(cols_to_use) == 1 && is.null(names(custom_colors))) {
    custom_colors <- stats::setNames(custom_colors, cols_to_use)
  } else if (is.list(custom_colors) && length(custom_colors) > 0 &&
    (is.null(names(custom_colors)) || any(names(custom_colors) == ""))) {
    cli::cli_abort("{.arg sample_colors} must be a named list when you supply colors for multiple annotations.")
  }

  for (column in cols_to_use) {
    cond_names <- unique(metadata[[column]])
    col_is_numeric <- ifelse(is.numeric(cond_names), TRUE, FALSE)
    cond_names[which(is.na(cond_names))] <- "NA"
    annot <- c()

    if (!is.null(custom_colors) && column %in% names(custom_colors)) {
      annot <- custom_colors[[column]]
      check_named_colors(annot, metadata[[column]], "sample_colors")
      annot_list[[column]] <- annot
      next
    }

    if (any(col_is_numeric, length(cond_names) > 9)) {
      palette <- grDevices::hcl.colors(length(cond_names), "Light Grays")
    } else {
      palette <- rep(default_palette, length.out = length(cond_names))
    }

    for (cond in cond_names) {
      annot[cond] <- palette[1]
      palette <- palette[-1]

    }
    annot_list[[column]] <- annot
  }
  return(annot_list)
}
