#' Regroup samples by condition and calculate mean for each condition
#'
#' @description An helper function that will give the mean by condition of the cscore or count, for each position. The standard deviation is also given.
#' this can be used to create a boxplot with ggplot.
#'
#' @param ribo a SummarizedExperiment object
#' @param metadata_condition name or index of the column __in metadata__ containing the condition
#' @param value name or index of the column containing the values on which mean by condition is calculated. (Either "cscore" or "count")
#' @param only_annotated use annotation site name instead of default
#'
#' @importFrom dplyr %>%
#' @importFrom rlang sym
#' @return a dataframe with the mean for each condition for a selected value
#' @export
#' @md
#' @examples
#' data("ribo_toy")
#' # mean_df <- mean_samples_by_condition(ribo_toy,"count","condition")
mean_samples_by_condition <- function(ribo, value, metadata_condition,
                                      only_annotated = FALSE) {
  val <- NULL
  check_is_se(ribo)
  check_type(value, "character", "value", length = 1)
  check_type(metadata_condition, "character", "metadata_condition", length = 1)
  check_type(only_annotated, "logical", "only_annotated", length = 1)
  check_metadata(ribo, metadata_condition)

  # Extract data in wide format
  # Returns [position_col, sample1, sample2...]
  df_wide <- extract_data(ribo, value, position_to_rownames = FALSE, only_annotated = only_annotated)

  # Identify position column name (first column)
  pos_col <- colnames(df_wide)[1]

  # Pivot longer
  # We assume all other columns are samples
  df_long <- tidyr::pivot_longer(df_wide,
    cols = -dplyr::all_of(pos_col),
    names_to = "samplename",
    values_to = "val"
  )

  # Add metadata
  # using colData
  cd <- as.data.frame(SummarizedExperiment::colData(ribo))
  cd$samplename <- rownames(cd) # ensure we have samplename to join

  # Join
  df_merged <- dplyr::left_join(df_long, cd, by = "samplename")

  # Summarize
  # Group by position and condition
  # Note: pos_col holds the site/position name. The output expects "site" as grouping var according to original code?
  # Original code: group_by(site, condition). "site" was likely the position column or explicitly 'site'.
  # I'll rename the position column to 'site' for consistency with downstream expectations if needed,
  # or just use pos_col. Use `site` if we want to match original output.

  grp_sym <- rlang::sym(pos_col)
  cond_sym <- rlang::sym(metadata_condition)

  ribo_condition <- df_merged %>%
    dplyr::group_by(!!grp_sym, !!cond_sym) %>%
    dplyr::summarise(mean = mean(val, na.rm = TRUE), sd = stats::sd(val, na.rm = TRUE))

  # Rename grouping column to 'site' if it is not
  if (pos_col != "site") {
    ribo_condition <- dplyr::rename(ribo_condition, site = !!grp_sym)
  }

  return(ribo_condition)
}
