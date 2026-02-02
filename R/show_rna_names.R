#' Display RNA names
#'
#' @param ribo a SummarizedExperiment object
#'
#' @return
#' A vector with actual RNA names
#' @export
#'
#' @examples
#' data("ribo_toy")
#' show_rna_names(ribo = ribo_toy)
show_rna_names <- function(ribo = NULL) {
  check_is_se(ribo)
  RNA_names <- SummarizedExperiment::metadata(ribo)$rna_names[[2]]
  return(RNA_names)
}
