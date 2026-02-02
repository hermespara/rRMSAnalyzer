#' Transform SummarizedExperiment data to a "ggplot-friendly" data frame.
#'
#' Turn a SummarizedExperiment into a dataframe that should be ready to be used in ggplot.
#' Metadata columns can be added if extra informations are needed.
#'
#' @param ribo A SummarizedExperiment object.
#' @param metadata_col Metadata columns to add. Must be in the object's metadata.
#' @param only_annotated Keep only sites that have been annotated.
#' @return A ggplot-friendly dataframe with the following columns : 1) site; 2) sample; 3) cscore; ...) any metadata added
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # format_to_plot(ribo_toy)
format_to_plot <- function(ribo, metadata_col = NULL, only_annotated = FALSE) {
  check_is_se(ribo)
  check_type(only_annotated, "logical", "only_annotated", length = 1)

  if (!is.null(metadata_col) && !is.character(metadata_col) && !is.numeric(metadata_col)) {
    cli::cli_abort("{.arg metadata_col} must be character or numeric.")
  }
  site <- NULL # NSE fix
  # First let's extract the data from ribo
  values_column <- "cscore"
  df.matrix <- extract_data(ribo, values_column,
    position_to_rownames = TRUE,
    only_annotated = only_annotated
  )
  # Add a new columns that correspond to rownames
  df.matrix$site <- rownames(df.matrix)
  # Transform the data with tidyr
  df.tranform <- tidyr::gather(df.matrix, "sample", "cscore", -site)

  if (is.null(metadata_col)) {
    # if no metadata, return the transformed data frame
    return(df.tranform)
  } else {
    # Check if the metadata_columns are numeric.
    # If true then get the column names of the selected columns

    cd <- as.data.frame(SummarizedExperiment::colData(ribo))

    if (is.numeric(metadata_col)) {
      metadata_columns <- colnames(cd[metadata_col])
    }
    # Get the metadata
    # Ensure we ask for existing columns
    cols_to_fetch <- c(metadata_col)

    # We need samplename to join. Assuming rownames of colData are samplenames.
    cd$samplename <- rownames(cd)

    # Select cols
    df.meta <- cd[, c("samplename", metadata_col), drop = FALSE]

    # if metadata is empty, return an error
    if (nrow(df.meta) == 0) {
      cli::cli_abort("No metadata found")
    }
    # Merge the metadata with the transformed data frame
    df.tranform <- merge(df.tranform, df.meta, by.x = "sample", by.y = "samplename")
    return(df.tranform)
  }
}
