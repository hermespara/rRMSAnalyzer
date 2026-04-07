# Generate the col argument for ComplexHeatmap::HeatmapAnnotation
generate_palette <- function(metadata, cols_to_use, custom_colors = NULL) {
  annot_list <- list()
  palettes_template <- list(
    c('#8dd3c7','#00bfff','#bebada','#f14292','#fdb462','#80b1d3','#b1de69','#fccde5','#d9d9d9'),
    c('#b2df8a','#ea1a8c','#a6cee3','#ffEf00','#1f78b4','#33a02c','#fddf6f','#cab2d6','#191970'),
    c('#fb9a99','#FFAA00','#AAFF00','#00FF00','#00FFAA','#00AAFF','#0000FF','#AA00FF','#FF00AA'))

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

  palettes <- palettes_template
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
      palette <- palettes[[1]]

    }

    for (cond in cond_names) {
      annot[cond] <- palette[1]
      palette <- palette[-1]

    }
    annot_list[[column]] <- annot
    palettes <- palettes[-1]
    if (length(palettes) == 0) palettes <- palettes_template
  }
  return(annot_list)
}
