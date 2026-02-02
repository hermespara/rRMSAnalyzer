#' Update count values inside a SummarizedExperiment with a matrix of values
#'
#' @param ribo a SummarizedExperiment object
#' @param matrix_new a position x sample matrix containing the new values
#'
#' @return a SummarizedExperiment with updated values
#' @keywords internal
#'
.update_ribo_count_with_matrix <- function(ribo, matrix_new) {
  check_is_se(ribo)
  check_type(matrix_new, "matrix", "matrix_new")

  # Verify dimensions
  if (!all(dim(matrix_new) == dim(ribo))) {
    cli::cli_abort("Dimensions of new matrix do not match SummarizedExperiment dimensions.")
  }

  # Ensure column names match (sanity check)
  if (!all(colnames(matrix_new) == colnames(ribo))) {
    warning("Column names of new matrix do not match SummarizedExperiment column names.")
  }

  SummarizedExperiment::assay(ribo, "counts", withDimnames = FALSE) <- matrix_new
  return(ribo)
}
