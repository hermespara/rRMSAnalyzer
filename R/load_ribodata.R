#' Load csv files from GenomeCov and their associated metadata. Create a SummarizedExperiment.
#' @inheritParams create_se
#' @inheritParams compute_cscore
#' @return a SummarizedExperiment object
#' @export
#'
#' @description Import your count CSV files and the metadata to create a SummarizedExperiment object.
#' The SummarizedExperiment object is used by rRMSAnalyzer package for all analyses.
#'
#' __This function serves as the entrypoint of rRMSAnalyzer.__
#'
#' @md
#' @details
#' load_ribodata is a wrapper of \code{\link{create_se}} and \code{\link{compute_cscore}}.
#' @seealso create_se
#'
load_ribodata <- function(count_path,
                          metadata = NULL,
                          count_sep = "\t",
                          metadata_sep = ",",
                          count_header = FALSE,
                          count_value = 3,
                          count_rnaid = 1,
                          count_pos = 2,
                          metadata_key = "filename",
                          metadata_id = NULL,
                          flanking = 6,
                          method = "median",
                          ncores = 1) {
  ribo <- create_se(
    count_path,
    metadata,
    count_sep,
    metadata_sep,
    count_header,
    count_value,
    count_rnaid,
    count_pos,
    metadata_key,
    metadata_id
  )

  check_type(flanking, "numeric", "flanking", length = 1)
  check_type(method, "character", "method", length = 1)
  check_type(ncores, "numeric", "ncores", length = 1)

  ribo <- compute_cscore(ribo, flanking, method, ncores)
  cli::cli_alert_success("Your data has been successfully loaded!")
  return(ribo)
}
