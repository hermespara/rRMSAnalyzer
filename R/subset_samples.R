#' Keep only selected samples from a SummarizedExperiment Object.
#'
#' @param ribo A SummarizedExperiment object.
#' @param samples_to_keep A vector containing the names of samples to keep.
#'
#' @return A SummarizedExperiment with only selected samples.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # ribo_test <- keep_samples(ribo = ribo_toy, samples_to_keep = c("S1", "S2"))
#'
keep_samples <- function(ribo = NULL, samples_to_keep = NULL) {
  check_is_se(ribo)
  check_type(samples_to_keep, "character", "samples_to_keep")
  check_sample(ribo, samples_to_keep)
  ribo <- ribo[, samples_to_keep]

  return(ribo)
}


#' Remove samples from a SummarizedExperiment object
#'
#' @param ribo A SummarizedExperiment object.
#' @param samples_to_delete A vector containing the names of samples to remove.
#'
#' @return A SummarizedExperiment without the removed samples.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # ribo_test <- remove_samples(ribo = ribo_toy, samples_to_delete = c("S1", "S2"))
#'
remove_samples <- function(ribo = NULL, samples_to_delete = NULL) {
  check_is_se(ribo)
  check_type(samples_to_delete, "character", "samples_to_delete")
  check_sample(ribo, samples_to_delete)

  # Ensure samples exist logic is in check_sample, but for removal we want to avoid error if some don't exist?
  # The original used match.
  # Let's keep strict check if check_sample provides it.

  idx <- match(samples_to_delete, colnames(ribo))
  if (any(is.na(idx))) {
    # If check_sample didn't catch it
    cli::cli_abort("Some samples to delete were not found in the object.")
  }
  ribo <- ribo[, -idx]

  return(ribo)
}
