#' Plot counts for a given position on a given RNA
#'
#'
#' @param ribo A SummarizedExperiment object
#' @param rna Name of RNA where the position is located
#' @param pos Position on RNA on which the view will be centered
#' @param samples Samples to display. "all" will display all samples.
#' @param flanking Number of sites to display on the left/right of the selected position.
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # ribo_toy <- rename_rna(ribo = ribo_toy)
#' # plot_counts_env(ribo = ribo_toy, rna = "5.8S", pos = 15)
plot_counts_env <- function(ribo = NULL, rna = NULL, pos = NULL, samples = "all", flanking = 6) {
  new_position <- count <- NULL

  check_is_se(ribo)
  check_type(rna, "character", "rna", length = 1)
  check_type(pos, "numeric", "pos", length = 1)
  check_type(flanking, "numeric", "flanking", length = 1)

  # check if rna is ok
  if (!(rna %in% ribo@metadata$rna_names[["current_name"]])) {
    cli::cli_abort("The RNA names given do not exist in the object")
  }

  rd <- SummarizedExperiment::rowData(ribo)



  # Calculate RNA length from rowData
  rna_len <- sum(rd$rna == rna)
  if (rna_len < pos) {
    stop(paste(pos, " is higher than the lenth of", rna))
  }

  # check for sample
  if (samples[1] == "all") {
    ribo <- ribo
  } else if (all(samples %in% colnames(ribo))) {
    ribo <- keep_samples(ribo = ribo, samples_to_keep = samples)
  } else {
    cli::cli_abort("Samples name should be from {toString(colnames(ribo))}")
  }

  # check for flanking
  # if(!(is.numeric(flanking))) {stop("flanking should be a number")}


  # format position of interest according to named_position of SummarizedExperiment
  pos_of_interest <- paste(rna, formatC(pos, width = 4, flag = "0"), sep = "_")

  # Extract count data from SummarizedExperiment
  # Returns dataframe with 'named_position' if position_to_rownames=F (default?)
  # extract_data defaults: col="cscore", pos_to_rn=F, only_annot=F.
  # We want "counts".

  count_data <- extract_data(ribo = ribo, col = "counts")

  # This returns df with [position, sample1, sample2...] where position is first column.
  # Assumption: check extract_data impl. It binds 'data' list.
  # In refactored extract_data: returns dataframe.
  # If position_to_rownames=FALSE, named_position column is present?
  # Need to verify extract_data return format.
  # Refactored extract_data:
  # matrix_all <- assay(ribo, col)
  # if !pos_to_rn: matrix_all <- cbind(named_position = rownames(matrix_all), as.data.frame(matrix_all))

  # extract the information
  which_pos <- which(count_data$named_position == pos_of_interest)

  if (length(which_pos) == 0) {
    stop(paste("Position", pos_of_interest, "not found in data."))
  }

  # positions around the position of interest
  # We need to handle boundary conditions? Original code didn't seem to explicitly?
  # It just slices.

  # Ensure indices are valid
  indices <- c((which_pos - flanking):which_pos, (which_pos + 1):(which_pos + flanking))
  indices <- indices[indices > 0 & indices <= nrow(count_data)]

  which_positions <- indices

  count_data <- count_data[which_positions, ]

  # new_position logic: this creates sequential x-axis.
  # If we are near boundaries, the count_data might be smaller than 2*flanking+1.
  # We should adjust new_position accordingly.

  # Re-calculate new_position vector based on captured rows
  # rowData has rnapos. We can use that?
  # Original: count_data$new_position <- c(...)
  # This assumes contiguous block centered on pos.

  # We can just map rnapos from rowData?
  # rd[which_positions, "rnapos"]

  count_data$new_position <- rd$rnapos[which_positions]

  # count_transform <- tidyr::gather(count_data[,-1], "samples", "count", -new_position)
  # count_data cols: named_position, S1, S2... new_position.
  # exclude named_position (col 1).

  count_transform <- tidyr::gather(count_data[, -1], "samples", "count", -new_position)


  # check if there are other modifications in the window
  # rd[which_positions, "site"]

  site_col <- rd$site[which_positions] # could be NA
  # Indices relative to the window
  # We want x-coordinates (new_position) where site is not NA.

  other_mod_indices <- which(!is.na(site_col))
  other_mod_pos <- count_data$new_position[other_mod_indices]


  # ggplot

  # Common plot logic...
  # Just copied and adapted:

  if (samples[1] == "all") {
    plot_to_return <- ggplot(data = count_transform) +
      geom_boxplot(aes(x = new_position, y = log10(count), group = new_position)) +
      # Highlight other mods
      {
        if (length(other_mod_pos) > 0) {
          geom_boxplot(
            data = count_transform[which(count_transform$new_position %in% other_mod_pos), ],
            aes(x = new_position, y = log10(count), group = new_position),
            fill = "#56B4E9", # Palette Blue
            width = 0.8
          )
        }
      } +
      # Highlight center
      geom_boxplot(
        data = count_transform[which(count_transform$new_position == pos), ],
        aes(x = new_position, y = log10(count), group = new_position),
        fill = "#009E73", # Palette Green
        width = 0.8
      ) +
      theme_rRMSAnalyzer() +
      labs(
        title = paste("Count profile for", ncol(ribo), "samples"),
        subtitle = paste("RNA:", rna),
        y = "log10(count)",
        x = "Position"
      ) +
      scale_x_continuous(
        labels = min(count_transform$new_position):max(count_transform$new_position),
        breaks = min(count_transform$new_position):max(count_transform$new_position)
      ) +
      {
        if (min(count_transform$count) < 100) {
          geom_hline(
            yintercept = 2,
            linewidth = 1,
            linetype = "dashed",
            color = "#D55E00"
          )
        }
      } + # Palette Vermilion
      {
        if (min(count_transform$count) < 100) {
          annotate("text",
            label = "Coverage limit",
            x = pos - flanking,
            y = 2 / 1.02,
            color = "#D55E00"
          )
        }
      } +
      geom_hline(
        yintercept = stats::median(log10(count_transform$count), na.rm = TRUE),
        linewidth = 1,
        linetype = "dashed",
        color = "#CC79A7"
      ) + # Palette Reddish Purple
      annotate("text",
        label = "Counts median",
        x = pos - flanking,
        y = stats::median(log10(count_transform$count) / 0.985, na.rm = TRUE),
        color = "#CC79A7"
      )
  } else {
    # samples specific
    plot_to_return <- ggplot(data = count_transform, aes(x = new_position, y = log10(count), group = samples)) +
      geom_point(size = 3) +
      geom_line(aes(col = samples), linewidth = 1.2) +
      theme_rRMSAnalyzer() +
      scale_color_rRMSAnalyzer() +
      labs(
        title = paste("Count profile for", length(samples), "samples"),
        subtitle = paste("RNA:", rna),
        y = "log10(count)",
        x = "position"
      ) +
      scale_x_continuous(
        labels = min(count_transform$new_position):max(count_transform$new_position),
        breaks = min(count_transform$new_position):max(count_transform$new_position)
      ) +
      {
        if (min(count_transform$count) < 100) {
          geom_hline(
            yintercept = 2,
            linewidth = 1,
            linetype = "dashed",
            color = "#D55E00"
          )
        }
      } +
      {
        if (min(count_transform$count) < 100) {
          annotate("text",
            label = "Coverage limit",
            x = pos - flanking + 1,
            y = 2 / 1.02,
            color = "#D55E00"
          )
        }
      } +
      geom_hline(
        yintercept = stats::median(log10(count_transform$count), na.rm = TRUE),
        linewidth = 1,
        linetype = "dashed",
        color = "#CC79A7"
      ) +
      annotate("text",
        label = "Counts median",
        x = pos - flanking + 1,
        y = stats::median(log10(count_transform$count) / 0.985, na.rm = TRUE),
        color = "#CC79A7"
      ) +
      geom_vline(
        xintercept = other_mod_pos,
        linewidth = 1,
        linetype = "dashed",
        color = "#56B4E9"
      ) +
      geom_vline(
        xintercept = pos,
        linewidth = 1,
        linetype = "dashed",
        color = "#009E73"
      )
  }

  return(plot_to_return)
}


# rajoute une ligne à 100
