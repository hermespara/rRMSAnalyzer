#' Correct batch effect with ComBat-seq method
#'
#' @description
#' Batch effect of RiboMethSeq data can be adjusted using the ComBat-seq method. adjust_bias is a wrapper to perform ComBat-seq adjustment.
#'
#' It will return a new SummarizedExperiment with adjusted read end count values and C-scores automatically recomputed with the same setup parameters.
#'
#'
#' @param ribo a SummarizedExperiment object.
#' @param batch Name of the column in metadata that contains the batch number.
#' @param ncores Number of cores to use in case of multithreading.
#' @param ... Parameters to pass to sva's \code{\link[sva]{ComBat_seq}} function.
#' @return SummarizedExperiment with adjusted read end count values and automatically recomputed C-scores.
#' @export
#'
#' @details
#' You must have a column with the batch number for each sample in your object's metadata (colData).
#'
#' @references Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
#'
#'
#' @examples
#' data("ribo_toy")
#' # ribo_toy_adjusted <- adjust_bias(ribo_toy,'run')
#'
adjust_bias <- function(ribo, batch, ncores = 1, ...) {
  check_is_se(ribo)
  check_type(batch, "character", "batch", length = 1)
  check_type(ncores, "numeric", "ncores", length = 1)
  check_metadata(ribo, batch)
  matrix_ribo <- extract_data(ribo, "counts", position_to_rownames = TRUE)

  # reorganize column according to metadata and convert DF to
  # matrix (otherwise, ComBat_seq won't work)
  # colData ensures alignment usually, but ComBat_seq might need strictly matrix
  col_data <- SummarizedExperiment::colData(ribo)

  # ComBat_seq expects matrix
  matrix_ribo <- as.matrix(matrix_ribo)

  adjusted_matrix <- sva::ComBat_seq(matrix_ribo,
    batch = col_data[[batch]],
    ...
  )
  ribo_updated <- .update_ribo_count_with_matrix(ribo, adjusted_matrix)

  if (isTRUE(ribo_updated@metadata$has_cscore)) {
    message(
      "Recomputing c-score with the following parameters :",
      "\n- C-score method : ", ribo_updated@metadata$cscore_method,
      "\n- Flanking window : ", ribo_updated@metadata$cscore_window,
      "\n"
    )
    ribo_updated <- compute_cscore(
      ribo_updated, ribo_updated@metadata$cscore_window,
      ribo_updated@metadata$cscore_method, ncores
    )
  }

  ribo_updated@metadata$combatSeq_count <- TRUE
  ribo_updated@metadata$col_used_combatSeq <- batch
  return(ribo_updated)
}
