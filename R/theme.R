#' Publication-ready theme for scientific journals
#'
#' A refined theme optimized for publication-quality figures. It provides
#' consistent typography, spacing, and panel styling across ggplot-based
#' visualizations in the package.
#'
#' @param base_size Base font size. Default is 11.
#' @param base_family Base font family. Default is \code{"sans"}.
#' @param base_line_size Base size for line elements. Default is 0.5.
#' @param base_rect_size Base size for rect elements. Default is 0.5.
#' @param aspect_ratio Optional aspect ratio (height/width). Use \code{NULL} or
#'   \code{"auto"} to leave the aspect unconstrained.
#' @return A ggplot2 theme object.
#' @import ggplot2
#' @export
theme_rRMSAnalyzer <- function(base_size = 11,
                               base_family = "sans",
                               base_line_size = 0.5,
                               base_rect_size = 0.5,
                               aspect_ratio = NULL) {
  if (identical(aspect_ratio, "auto")) {
    aspect_ratio <- NULL
  }

  if (!is.null(aspect_ratio)) {
    if (!is.numeric(aspect_ratio) || length(aspect_ratio) != 1 || is.na(aspect_ratio)) {
      cli::cli_abort("{.arg aspect_ratio} must be a single numeric value, {.val NULL}, or {.val 'auto'}.")
    }
  }

  if (base_family == "sans") {
    if (.Platform$OS.type == "windows") {
      base_family <- "Arial"
    } else if (Sys.info()[["sysname"]] == "Darwin") {
      base_family <- "Helvetica"
    }
  }

  half_line <- base_size / 2

  theme_custom <- ggplot2::theme_bw(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.15),
        hjust = 0,
        margin = ggplot2::margin(b = half_line)
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(0.95),
        color = "grey30",
        hjust = 0,
        margin = ggplot2::margin(b = half_line)
      ),
      plot.caption = ggplot2::element_text(
        size = ggplot2::rel(0.85),
        color = "grey50",
        hjust = 1,
        margin = ggplot2::margin(t = half_line)
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.2),
        hjust = 0,
        vjust = 1
      ),
      plot.tag.position = c(0, 1),
      axis.title.x = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.0),
        margin = ggplot2::margin(t = half_line * 0.8)
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.0),
        margin = ggplot2::margin(r = half_line * 0.8),
        angle = 90
      ),
      axis.text.x = ggplot2::element_text(
        color = "black",
        size = ggplot2::rel(0.95),
        margin = ggplot2::margin(t = 0.2 * base_size)
      ),
      axis.text.y = ggplot2::element_text(
        color = "black",
        size = ggplot2::rel(0.95),
        margin = ggplot2::margin(r = 0.2 * base_size)
      ),
      axis.ticks = ggplot2::element_line(color = "grey20", linewidth = 0.5),
      axis.ticks.length = ggplot2::unit(0.15, "cm"),
      axis.line = ggplot2::element_blank(),
      legend.position = "right",
      legend.justification = c(0, 1),
      legend.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(0.95)
      ),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      legend.spacing = ggplot2::unit(0.4, "cm"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(0.4, "cm"),
      strip.background = ggplot2::element_rect(
        fill = "grey90",
        color = "grey20",
        linewidth = 0.5
      ),
      strip.text.x = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(0.95),
        margin = ggplot2::margin(0.15 * base_size, 0, 0.15 * base_size, 0)
      ),
      strip.text.y = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(0.95),
        angle = -90,
        margin = ggplot2::margin(0, 0.15 * base_size, 0, 0.15 * base_size)
      ),
      strip.placement = "outside",
      strip.switch.pad.grid = ggplot2::unit(0.1, "cm"),
      strip.switch.pad.wrap = ggplot2::unit(0.1, "cm"),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(
        fill = NA,
        color = "grey20",
        linewidth = 0.7
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.4
      ),
      panel.grid.major.y = ggplot2::element_line(
        color = "grey92",
        linewidth = 0.4
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(half_line, "pt"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(half_line, half_line, half_line, half_line)
    )

  if (!is.null(aspect_ratio)) {
    theme_custom <- theme_custom + ggplot2::theme(aspect.ratio = aspect_ratio)
  }

  theme_custom
}

# Enhanced color palette based on publication-friendly, colorblind-safe colors.
.rRMSAnalyzer_palette <- function() {
  c(
    "#0077BB",
    "#CC3311",
    "#009988",
    "#EE7733",
    "#33BBEE",
    "#EE3377",
    "#BBBBBB",
    "#000000",
    "#117733",
    "#882255",
    "#44AA99",
    "#88CCEE",
    "#DDCC77",
    "#AA4499",
    "#DDDDDD",
    "#332288"
  )
}

# Grayscale palette for journal requirements.
.rRMSAnalyzer_grayscale <- function() {
  c(
    "#000000", "#404040", "#808080", "#BFBFBF",
    "#E0E0E0", "#F0F0F0", "#555555", "#999999"
  )
}

# Scientific diverging palette for continuous values.
.rRMSAnalyzer_diverging_palette <- function(palette = "blue_red", direction = 1) {
  colors <- switch(palette,
    blue_red = c(
      "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
      "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
    ),
    blue_white_red = c(
      "#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020"
    ),
    NULL
  )

  if (!is.null(colors) && identical(direction, -1)) {
    colors <- rev(colors)
  }

  colors
}

#' Custom color scale for rRMSAnalyzer
#'
#' Colorblind-friendly palette optimized for scientific publications.
#'
#' @param ... Arguments passed to \code{\link[ggplot2]{scale_color_manual}}.
#' @param values Optional character vector of colors. If \code{NULL}, uses the
#'   rRMSAnalyzer palette.
#' @param grayscale Logical. If \code{TRUE}, uses the grayscale palette.
#' @return A discrete color scale.
#' @import ggplot2
#' @export
scale_color_rRMSAnalyzer <- function(..., values = NULL, grayscale = FALSE) {
  if (is.null(values)) {
    values <- if (isTRUE(grayscale)) {
      .rRMSAnalyzer_grayscale()
    } else {
      .rRMSAnalyzer_palette()
    }
  }

  ggplot2::scale_color_manual(values = values, ...)
}

#' Custom fill scale for rRMSAnalyzer
#'
#' Colorblind-friendly palette optimized for scientific publications.
#'
#' @param ... Arguments passed to \code{\link[ggplot2]{scale_fill_manual}}.
#' @param values Optional character vector of colors. If \code{NULL}, uses the
#'   rRMSAnalyzer palette.
#' @param grayscale Logical. If \code{TRUE}, uses the grayscale palette.
#' @return A discrete fill scale.
#' @import ggplot2
#' @export
scale_fill_rRMSAnalyzer <- function(..., values = NULL, grayscale = FALSE) {
  if (is.null(values)) {
    values <- if (isTRUE(grayscale)) {
      .rRMSAnalyzer_grayscale()
    } else {
      .rRMSAnalyzer_palette()
    }
  }

  ggplot2::scale_fill_manual(values = values, ...)
}

#' Continuous color scale for heatmaps and gradients
#'
#' A perceptually uniform continuous scale suitable for scientific figures.
#'
#' @param ... Arguments passed to \code{\link[ggplot2]{scale_color_gradientn}}.
#' @param palette Character string specifying the palette. One of
#'   \code{"viridis"}, \code{"plasma"}, \code{"inferno"}, \code{"magma"},
#'   \code{"cividis"}, \code{"blue_red"}, or \code{"blue_white_red"}.
#' @param direction Sets the order of colors. \code{1} (default) or
#'   \code{-1} (reversed).
#' @return A continuous color scale.
#' @import ggplot2
#' @export
scale_color_rRMSAnalyzer_c <- function(..., palette = "viridis", direction = 1) {
  if (palette %in% c("viridis", "plasma", "inferno", "magma", "cividis")) {
    return(ggplot2::scale_color_viridis_c(option = palette, direction = direction, ...))
  }

  colors <- .rRMSAnalyzer_diverging_palette(palette = palette, direction = direction)
  if (is.null(colors)) {
    return(ggplot2::scale_color_viridis_c(direction = direction, ...))
  }

  ggplot2::scale_color_gradientn(colors = colors, ...)
}

#' Continuous fill scale for heatmaps and gradients
#'
#' A perceptually uniform continuous scale suitable for scientific figures.
#'
#' @param ... Arguments passed to \code{\link[ggplot2]{scale_fill_gradientn}}.
#' @param palette Character string specifying the palette. One of
#'   \code{"viridis"}, \code{"plasma"}, \code{"inferno"}, \code{"magma"},
#'   \code{"cividis"}, \code{"blue_red"}, or \code{"blue_white_red"}.
#' @param direction Sets the order of colors. \code{1} (default) or
#'   \code{-1} (reversed).
#' @return A continuous fill scale.
#' @import ggplot2
#' @export
scale_fill_rRMSAnalyzer_c <- function(..., palette = "viridis", direction = 1) {
  if (palette %in% c("viridis", "plasma", "inferno", "magma", "cividis")) {
    return(ggplot2::scale_fill_viridis_c(option = palette, direction = direction, ...))
  }

  colors <- .rRMSAnalyzer_diverging_palette(palette = palette, direction = direction)
  if (is.null(colors)) {
    return(ggplot2::scale_fill_viridis_c(direction = direction, ...))
  }

  ggplot2::scale_fill_gradientn(colors = colors, ...)
}

#' Save publication-ready figures
#'
#' Saves figures with journal-appropriate dimensions and resolution.
#'
#' @param filename File name to save, including the extension.
#' @param plot ggplot object. If \code{NULL}, saves the last plot.
#' @param width Width in mm. Common values are 89 (single column) and 183
#'   (double column).
#' @param height Height in mm.
#' @param dpi Resolution in dots per inch. Default is 300.
#' @param format Optional journal preset. One of \code{"nature_single"},
#'   \code{"nature_double"}, \code{"pnas_single"}, \code{"pnas_double"},
#'   \code{"science_single"}, \code{"science_double"}, or \code{"custom"}.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{ggsave}}.
#' @return Invisibly returns the saved file name.
#' @export
save_publication_figure <- function(filename,
                                    plot = ggplot2::last_plot(),
                                    width = NULL,
                                    height = NULL,
                                    dpi = 300,
                                    format = "custom",
                                    ...) {
  dimensions <- list(
    nature_single = list(width = 89, height = 89),
    nature_double = list(width = 183, height = 183),
    pnas_single = list(width = 87, height = 87),
    pnas_double = list(width = 178, height = 178),
    science_single = list(width = 90, height = 90),
    science_double = list(width = 183, height = 183)
  )

  if (format != "custom" && format %in% names(dimensions)) {
    if (is.null(width)) {
      width <- dimensions[[format]]$width
    }
    if (is.null(height)) {
      height <- dimensions[[format]]$height
    }
  }

  if (is.null(width)) {
    width <- 89
  }
  if (is.null(height)) {
    height <- 89
  }

  ext <- tools::file_ext(filename)
  device <- switch(ext,
    pdf = "pdf",
    png = "png",
    tiff = "tiff",
    tif = "tiff",
    eps = "eps",
    "pdf"
  )

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "mm",
    dpi = dpi,
    device = device,
    ...
  )

  message(sprintf(
    "Figure saved: %s (%d x %d mm at %d dpi)",
    filename, width, height, dpi
  ))

  invisible(filename)
}
