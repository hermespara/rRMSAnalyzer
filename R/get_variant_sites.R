#' Return the most/less variant sites of a dataframe or a SummarizedExperiment
#'
#' @param ribo a SummarizedExperiment object.
#' @param n Number of top sites to return.
#' @param type_of_variant Either "less" or "most", to select the top n less
#' variant sites or the top n most variant sites respectively.
#' @param only_annotated (SummarizedExperiment only) Check variability only among
#' annotated sites. Ignored when df is a dataframe.
#'
#' @return A dataframe with the n most/less variant sites
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' ribo_toy <- rename_rna(ribo_toy)
#' ribo_toy <- annotate_site(ribo_toy, human_methylated)
#' get_variant_sites(ribo = ribo_toy)
get_variant_sites <- function(ribo, n = 20, type_of_variant = "most",
                              only_annotated = TRUE) {
  check_is_se(ribo)
  check_type(n, "numeric", "n", length = 1)
  check_type(type_of_variant, "character", "type_of_variant", length = 1)
  check_in_set(tolower(type_of_variant), c("most", "less"), "type_of_variant")
  check_type(only_annotated, "logical", "only_annotated", length = 1)

  site <- NULL

  df <- extract_data(ribo,
    only_annotated = only_annotated,
    position_to_rownames = TRUE
  )

  df <- as.data.frame(df)

  if (tolower(type_of_variant) == "most") {
    most_variant <- TRUE
  } else {
    most_variant <- FALSE
  }

  var <- apply(df, 1, var)
  df <- df[order(var, decreasing = most_variant)[1:n], ]
  df["site"] <- rownames(df)
  df <- dplyr::relocate(df, site)
  return(df)
}
