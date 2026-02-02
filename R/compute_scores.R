# (internal) Compute a c-score for a vector of counts for a single RNA
# (internal) Compute a c-score for a vector of counts for a single RNA
.compute_vec_cscore <- function(count, flanking = 6, method) {
  # Call Rcpp function to get flanking (window) values
  flanking_values <- compute_window_score_cpp(as.numeric(count), as.numeric(flanking), as.character(method))

  if (method %in% c("median", "mean")) {
    scorec_raw <- 1 - count / flanking_values
  }
  return(pmax(scorec_raw, 0))
}

# (internal) Compute c-score for one sample (all RNAs)
.compute_sample_cscore_se <- function(sample_counts, rna_factor, flanking, method) {
  # split counts by RNA
  counts_by_rna <- split(sample_counts, rna_factor)

  # compute score for each RNA
  scores_by_rna <- lapply(counts_by_rna, .compute_vec_cscore, flanking, method)

  # unsplit to get back original order
  scores <- unsplit(scores_by_rna, rna_factor)
  return(scores)
}

#' Compute c-score for all samples in a SummarizedExperiment
#'
#' @description
#'
# The C-score corresponds to the 2'Ome level at a RNA position. The
# C-score represents a drop in the end read coverage at a given
# position compared to the environmental coverage, as described by
# @birkedal2014. The C-score can be of 0 (i.e., no RNA molecule is
# 2'Ome at the position of interest), of 1 (i.e., all the RNA
# molecules are 2'Ome at the position of interest) and of ]0:1[
# (i.e., a mix of un-methylated and methylated RNA molecules).
#'
#' In this package, the C-score is calculated for every position. As it can be
#' useful to find positions not yet identified as methylated.
#'
#' This function is called automatically when loading data or after adjusting
#' biases in end read count data with ComBat-seq. It can be still called
#' manually to modify how the C-score is currently computed.
#'
#' For each RNA, the first and last positions cannot be calculated if the local
#' coverage is shorter than the flanking argument. Their value will be NA instead.
#'
#' @references Birkedal, U., Christensen-Dalsgaard, M., Krogh, N., Sabarinathan,
#' R., Gorodkin, J. and Nielsen, H. (2015),
#' Profiling of Ribose Methylations in RNA by High-Throughput Sequencing.
#' Angew. Chem. Int. Ed., 54: 451-455. https://doi.org/10.1002/anie.201408362
#'
#'
#' @md
#' @param ribo A SummarizedExperiment object.
#' @param flanking Size of the local coverage.
#' @param method Computation method of the local coverage. Either 'median' or 'mean'.
#' @param ncores Number of cores to use in case of multithreading.
#' @return A SummarizedExperiment with c-score assay added.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # ribo_with_cscore_med <- compute_cscore(ribo_toy, ncores = 2)
#'
compute_cscore <- function(ribo = NULL, flanking = 6, method = "median",
                           ncores = 1) {
  check_is_se(ribo)
  check_type(flanking, "numeric", "flanking", length = 1)
  check_type(method, "character", "method", length = 1)
  check_in_set(method, c("median", "mean"), "method")
  check_type(ncores, "numeric", "ncores", length = 1)

  counts <- SummarizedExperiment::assay(ribo, "counts")
  rnas <- SummarizedExperiment::rowData(ribo)$rna

  # Use parallel processing or standard apply
  if (ncores > 1) {
    # Parallel apply over columns (samples)
    # We need to transpose the result of mclapply if iterating over cols?
    # Actually mclapply over columns of a matrix isn't direct.
    # We can use index.
    sample_indices <- seq_len(ncol(counts))
    cscores_list <- parallel::mclapply(sample_indices, function(i) {
      .compute_sample_cscore_se(counts[, i], rnas, flanking, method)
    }, mc.cores = ncores)
    cscores <- do.call(cbind, cscores_list)
    colnames(cscores) <- colnames(counts)
  } else {
    cscores <- apply(counts, 2, .compute_sample_cscore_se, rna_factor = rnas, flanking = flanking, method = method)
  }

  SummarizedExperiment::assay(ribo, "cscore") <- cscores
  ribo@metadata$cscore_window <- flanking
  ribo@metadata$cscore_method <- method
  ribo@metadata$has_cscore <- TRUE

  return(ribo)
}
