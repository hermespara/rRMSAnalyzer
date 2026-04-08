#' Remove one or more RNA among your SummarizedExperiment' samples
#'
#' @param ribo A SummarizedExperiment object.
#' @param rna_to_remove One RNA name or a vector of RNA names.
#'
#' @return A SummarizedExperiment with the RNA(s) removed.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # ribo_toy <- remove_rna(ribo_toy, "NR_023363.1_5S")
remove_rna <- function(ribo, rna_to_remove) {
  check_is_se(ribo)
  check_type(rna_to_remove, "character", "rna_to_remove")

  metadata_list <- S4Vectors::metadata(ribo)
  rna_names_df <- metadata_list$rna_names
  if (!all(rna_to_remove %in% rna_names_df[["current_name"]])) {
    cli::cli_abort("The RNA names given do not exist in the object")
  }

  # Drop unused levels for 'rna' factor in rowData is important?
  # For efficiency, we just subset. SE handles it.
  # But we might want to drop levels if we access it later.

  rows_to_keep <- !(SummarizedExperiment::rowData(ribo)$rna %in% rna_to_remove)
  ribo <- ribo[rows_to_keep, ]

  # Drop levels in rowData$rna
  # Note: rowData is a DataFrame. Modifying it:
  rd <- SummarizedExperiment::rowData(ribo)
  rd$rna <- droplevels(rd$rna)
  SummarizedExperiment::rowData(ribo) <- rd

  # Update internal metadata for rna names
  metadata_list$rna_names <- rna_names_df[!(rna_names_df[["current_name"]] %in% rna_to_remove), ]
  S4Vectors::metadata(ribo) <- metadata_list

  return(ribo)
}
