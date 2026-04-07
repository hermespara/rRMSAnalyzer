#' Unified theme for rRMSAnalyzer plots
#'
#' A publication-ready theme for ggplot2 plots, ensuring consistency across all
#' visualizations in the package. It features a clean, minimal design with
#' clear typography and distinctive lines.
#'
#' @param base_size Base font size. Default is 14.
#' @param base_family Base font family. Default is "sans".
#' @return A ggplot2 theme.
#' @import ggplot2
#' @export
theme_rRMSAnalyzer <- function(base_size = 14, base_family = "sans") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
            # Text
            plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0, margin = ggplot2::margin(b = 10)),
            plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.95), color = "grey30", hjust = 0, margin = ggplot2::margin(b = 10)),
            plot.caption = ggplot2::element_text(size = ggplot2::rel(0.8), color = "grey50", hjust = 1, margin = ggplot2::margin(t = 10)),
            axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0)),
            axis.text = ggplot2::element_text(color = "black", size = ggplot2::rel(0.9)),

            # Legend
            legend.position = "right",
            legend.title = ggplot2::element_text(face = "bold"),
            legend.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank(),

            # Panels and Grids
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.5),
            panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.5),
            panel.border = ggplot2::element_rect(fill = NA, color = "grey20", linewidth = 0.8),

            # Background
            plot.background = ggplot2::element_rect(fill = "white", color = NA),
            panel.background = ggplot2::element_rect(fill = "white", color = NA)
        )
}

#' Custom color palette for rRMSAnalyzer
#'
#' A colorblind-friendly palette for categorical variables.
#'
#' @param ... Arguments passed to \code{\link[ggplot2]{scale_color_manual}}.
#' @param values Optional named character vector of colors to use instead of
#' the default rRMSAnalyzer palette.
#' @return A discrete color scale.
#' @import ggplot2
#' @export
scale_color_rRMSAnalyzer <- function(..., values = NULL) {
    # Okabe-Ito palette (Colorblind friendly)
    palette <- values
    if (is.null(palette)) {
        palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    }
    ggplot2::scale_color_manual(values = palette, ...)
}

#' Custom fill palette for rRMSAnalyzer
#'
#' A colorblind-friendly palette for categorical variables.
#'
#' @param ... Arguments passed to \code{\link[ggplot2]{scale_fill_manual}}.
#' @param values Optional named character vector of colors to use instead of
#' the default rRMSAnalyzer palette.
#' @return A discrete fill scale.
#' @import ggplot2
#' @export
scale_fill_rRMSAnalyzer <- function(..., values = NULL) {
    palette <- values
    if (is.null(palette)) {
        palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
    }
    ggplot2::scale_fill_manual(values = palette, ...)
}
