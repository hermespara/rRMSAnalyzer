#' Aggregate results into a single matrix
#'
#' For a given column in data, this function will generate a dataframe with all samples.
#' Exports all positions (if only_annotated is false) or only annotated sites
#' (if only_annotated is true).
#'
#' @param ribo A SummarizedExperiment object.
#' @param col Assay name to extract (e.g., "cscore" or "counts").
#' @param position_to_rownames If true, position will be included as a rowname.
#' They will in a new column otherwise.
#' @param only_annotated If true, return a dataframe with only annotated sites.
#' Return all sites otherwise.
#' @return A position x samples dataframe. Position can be all or only
#' annotated sites. A site column can be added if position_to_rownames = FALSE.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # count_df <- extract_data(ribo_toy,'count')
extract_data <- function(ribo, col = "cscore",
                         position_to_rownames = FALSE, only_annotated = FALSE) {
  check_is_se(ribo)
  check_type(col, "character", "col", length = 1)
  check_type(position_to_rownames, "logical", "position_to_rownames", length = 1)
  check_type(only_annotated, "logical", "only_annotated", length = 1)

  col <- tolower(col)
  if (col == "count") col <- "counts" # Alias handling

  if (!(col %in% SummarizedExperiment::assayNames(ribo))) {
    cli::cli_abort(c("Name supplied to {.var col} is not an assay in the object",
      "i" = "Available assays: {.val {SummarizedExperiment::assayNames(ribo)}}",
      "x" = "{.val {col}} is not a valid assay"
    ))
  }

  # Extract the assay data
  mat <- SummarizedExperiment::assay(ribo, col)
  df <- as.data.frame(mat)

  # Get row details
  rd <- SummarizedExperiment::rowData(ribo)

  # Handle "only_annotated"
  if (only_annotated) {
    if (!"site" %in% names(rd)) {
      cli::cli_abort("No 'site' column in rowData. Can't filter by annotated sites.")
    }
    keep <- !is.na(rd$site) & rd$site != ""
    df <- df[keep, , drop = FALSE]
    rd <- rd[keep, , drop = FALSE]
    position_col_name <- "site"
    position_values <- rd$site
  } else {
    # Generate named_position if not present or use constructed one
    # Original logic used named_position which was likely "RNA_Position"
    # We can reconstruct it or check if it's in rowData
    if ("named_position" %in% names(rd)) {
      position_values <- rd$named_position
    } else {
      position_values <- paste(rd$rna, formatC(rd$rnapos, width = 4, flag = "0"), sep = "_")
    }
    position_col_name <- "named_position"
  }

  # Check for match (logic from original: match(position_list, matrix_all...))
  # In SE, strict alignment is guaranteed, so we don't need to join/match like before.

  if (position_to_rownames) {
    rownames(df) <- position_values
  } else {
    # Add position column at the beginning
    df <- cbind(setNames(data.frame(position_values), position_col_name), df)
  }

  return(df)
}
