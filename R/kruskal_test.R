#' Kruskal test on c-scores
#'
#' @param cscore_matrix Sites x Samples C-score matrix (output of extract_data()).
#' @param metadata Metadata associated to samples c-score matrix.
#' @param adjust_pvalues_method Method to adjust p-value, one of p.adjust.methods.
#' @param factor_column Metadata column used to group samples by.
#' @param statistical_test Statistical test used to compute p-values.
#' One of "kruskal" (wilcox is automatically used if there are 2 groups) or "t.test".
#'
#' @return A dataframe with three columns : 1) site : RNA site, 2) p.val : raw p-values, 3) adjusted p-values.
#' @keywords internal
#'
kruskal_test_on_cscores <- function(cscore_matrix = NULL,
                                    metadata = NULL,
                                    adjust_pvalues_method = "fdr",
                                    factor_column = NULL,
                                    statistical_test = "kruskal") {
  check_type(cscore_matrix, "data.frame", "cscore_matrix")
  check_type(metadata, "data.frame", "metadata")
  check_type(adjust_pvalues_method, "character", "adjust_pvalues_method", length = 1)
  check_type(factor_column, "character", "factor_column", length = 1)
  check_type(statistical_test, "character", "statistical_test", length = 1)
  check_in_set(tolower(statistical_test), c("kruskal", "t.test"), "statistical_test")

  if (!factor_column %in% colnames(metadata)) {
    cli::cli_abort("{.arg factor_column} {.val {factor_column}} is not present in metadata.")
  }

  cscore_matrix <- as.data.frame(cscore_matrix)
  cscore_matrix <- cscore_matrix[stats::complete.cases(cscore_matrix), match(metadata[, "samplename"], colnames(cscore_matrix))] # order column cscores as metadata
  n_groups <- length(unique(metadata[, factor_column]))
  test_to_perform <- tolower(statistical_test)

  if (test_to_perform == "kruskal") {
    test_to_perform <- ifelse(n_groups > 2, "kruskal", "wilcox")
  }
  if (test_to_perform == "t.test" && n_groups != 2) {
    cli::cli_abort("{.arg statistical_test} = {.val t.test} requires exactly 2 groups in {.arg factor_column}. Found {n_groups}.")
  }

  if (test_to_perform == "kruskal") {
    kruskal_test_pvalues <- apply(cscore_matrix, 1, function(x) {
      stats::kruskal.test(x ~ metadata[, factor_column])$p.value
    })
  }
  if (test_to_perform == "wilcox") {
    kruskal_test_pvalues <- apply(cscore_matrix, 1, function(x) {
      stats::wilcox.test(x ~ metadata[, factor_column], exact = TRUE)$p.value # test each row and get p.value
    })
  }
  if (test_to_perform == "t.test") {
    kruskal_test_pvalues <- apply(cscore_matrix, 1, function(x) {
      stats::t.test(x ~ metadata[, factor_column])$p.value
    })
  }

  df_kruskal_pvalues <- data.frame(
    site = names(kruskal_test_pvalues), p.val = kruskal_test_pvalues,
    p.adj = stats::p.adjust(kruskal_test_pvalues, method = adjust_pvalues_method)
  )

  # Return a data frame with both raw and adjusted p-values.
  # Example :
  #   site      p.val     p.adj
  # 18S_Am27  0.5085506 0.6929003
  # 18S_Am99  0.4942262 0.6929003
  # 18S_Um116 0.0120457 0.0883697
  # 18S_Um121 0.4803053 0.6888589
  # 18S_Am159 0.2537442 0.4687817
  return(df_kruskal_pvalues)
}



#' Return the range between conditions' mean
#'
#' @inheritParams kruskal_test_on_cscores
#'
#' @keywords internal
#' @return a dataframe with 1) annotated site name 2) range of mean between conditions for this site.
#'
get_range <- function(cscore_matrix = NULL, metadata = NULL, factor_column = NULL) {
  # TODO : rename this function ?

  check_type(cscore_matrix, "data.frame", "cscore_matrix")
  check_type(metadata, "data.frame", "metadata")
  check_type(factor_column, "character", "factor_column", length = 1)

  if (!factor_column %in% colnames(metadata)) {
    cli::cli_abort("{.arg factor_column} {.val {factor_column}} is not present in metadata.")
  }

  cscore_matrix <- t(cscore_matrix)
  cscore_matrix <- as.data.frame(cscore_matrix[
    match(metadata[, "samplename"], rownames(cscore_matrix)),
  ])

  # Aggregate each site, compute the mean according to factor_column
  df_mean_each_group <- stats::aggregate(
    cscore_matrix,
    list(metadata[, factor_column]), mean
  )

  # compute the difference between the max - min group
  df_min_max <- apply(df_mean_each_group[, -1], 2, function(x) {
    diff(range(x))
  })

  df_min_max <- data.frame(
    site = names(df_min_max),
    mean_max_min_difference = df_min_max
  ) # df summary


  # Example of return :
  #   site          mean_max_min_difference
  # 18S_Am27             0.003995770
  # 18S_Am99             0.004179315
  # 18S_Um116            0.039218474
  # 18S_Um121            0.005993061
  # 18S_Am159            0.007420624

  return(df_min_max)
}

#' (internal) Test which sites are the most variable between conditions
#'
#' This function outputs the same result as plot_diff_sites, but as a dataframe.
#' It cannot be called directly, use "only_df = TRUE" in plot_diff_sites instead.
#'
#' @inherit plot_diff_sites details
#'
#' @param ribo a SummarizedExperiment
#' @param adjust_pvalues_method p-value adjustment method (default : "fdr")
#' @param factor_column Metadata column used to group samples by.
#' @param statistical_test Statistical test used to compute p-values.
#' One of "kruskal" (wilcox is automatically used if there are 2 groups) or "t.test".
#'
#' @return A dataframe where each site is a row and with the following column :
#'  1) site : name of the site; 2) p.val : raw p-value; 3) p.adj : adjusted p-value 4) mean_max_min_difference : range between conditions of mean c-score
#' @keywords internal
wrapper_kruskal_test <- function(ribo = NULL,
                                 adjust_pvalues_method = "fdr",
                                 factor_column = NULL,
                                 statistical_test = "kruskal") {
  check_is_se(ribo)
  check_type(adjust_pvalues_method, "character", "adjust_pvalues_method", length = 1)
  check_type(factor_column, "character", "factor_column", length = 1)
  check_type(statistical_test, "character", "statistical_test", length = 1)
  check_in_set(tolower(statistical_test), c("kruskal", "t.test"), "statistical_test")
  check_metadata(ribo, factor_column)

  metadata <- as.data.frame(SummarizedExperiment::colData(ribo))
  # Ensure samplename exists
  if (!"samplename" %in% colnames(metadata)) {
    metadata$samplename <- rownames(metadata)
  }

  cscore_matrix <- extract_data(ribo, only_annotated = TRUE, position_to_rownames = TRUE)
  if (nrow(cscore_matrix) == 0) {
    if (nrow(cscore_matrix) == 0) {
      cli::cli_abort("Please annotate your object before using plot_diff_sites")
    }
  }
  df_pval <- kruskal_test_on_cscores(
    cscore_matrix = cscore_matrix,
    metadata = metadata,
    adjust_pvalues_method = adjust_pvalues_method,
    factor_column = factor_column,
    statistical_test = statistical_test
  )

  df_min_max <- get_range(
    cscore_matrix = cscore_matrix,
    metadata = metadata,
    factor_column = factor_column
  )

  df_final <- merge(df_pval, df_min_max, by = "site")

  # Example of df_final :
  #     site     p.val     p.adj         mean_max_min_difference
  # 18S_Am1031 0.2306932 0.4335441             0.004622258
  # 18S_Am1383 0.4408513 0.6582574             0.003166407
  # 18S_Am159  0.2537442 0.4687817             0.007420624
  # 18S_Am166  0.2158151 0.4277063             0.008474824
  # 18S_Am1678 0.2220702 0.4322437             0.015554299

  return(df_final)
}
