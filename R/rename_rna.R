#' Rename RNAs in a SummarizedExperiment object
#'
#' @param ribo A SummarizedExperiment object.
#' @param new_names Vector of new RNA names (by RNA size order).
#' @return a SummarizedExperiment with updated RNA names.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # ribo_toy <- rename_rna(ribo_toy ,c("5S","5.8S","18S","28S"))
#'
rename_rna <- function(ribo, new_names = c("5S", "5.8S", "18S", "28S")) {
  check_is_se(ribo)
  check_type(new_names, "character", "new_names")

  metadata_list <- S4Vectors::metadata(ribo)
  rna_names_df <- metadata_list$rna_names

  if (nrow(rna_names_df) != length(new_names)) {
    cli::cli_abort("Different numbers of RNA names in your object ({nrow(rna_names_df)}) and the list given ({length(new_names)}).")
  }

  # Create mapping
  # The rna_names_df has rows corresponding to distinct RNAs.
  # We assume the order of new_names corresponds to the order of rows in rna_names_df?
  # The doc says "by RNA size order".
  # Original code: rna_names[3] <- new_names. It assumed rna_names rows were ordered.

  # Update rowData
  rd <- SummarizedExperiment::rowData(ribo)

  # Map current names to new names
  current_names <- as.character(rd$rna)

  # Map:
  # rna_names_df$current_name[i] -> new_names[i]

  # We can use factor levels replacement if we align them
  mapping <- stats::setNames(new_names, rna_names_df$current_name)

  new_rna_col <- mapping[current_names]

  # Update rd$rna
  rd$rna <- factor(new_rna_col, levels = unique(new_names)) # Maintain order?

  SummarizedExperiment::rowData(ribo) <- rd

  # Update metadata
  rna_names_df$current_name <- new_names
  metadata_list$rna_names <- rna_names_df
  S4Vectors::metadata(ribo) <- metadata_list

  return(ribo)
}
