#' Annotate sites according to a dataframe containing annotation
#'
#' Annotate sites of interest by giving them a custom name.
#' Some analyses and plots can be applied on these specific sites when only_annotated parameter
#' is available.
#'
#' @details
#' This function will fill the 'site' column in your object's rowData with a nomenclature given in annot.
#'
#' @param ribo a SummarizedExperiment object to annotate, see :
#' \code{\link{load_ribodata}}
#' @param annot The dataframe containing annotations.
#' @param anno_rna Name or index of the column in annot containing RNAs' name.
#' @param anno_pos Name or index of the column in annot containing site position inside RNA.
#' @param anno_value Name or index of the column in annot containing nomenclature to apply.
#'
#' @return An annotated SummarizedExperiment (the site column in rowData should be filled with names from annot_value for known position).
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' # ribo_toy <- rename_rna(ribo_toy ,c('5S','5.8S','18S','28S'))
#' # ribo_toy <- annotate_site(ribo_toy,human_methylated,anno_value ='Nomenclature')
#'
annotate_site <- function(ribo, annot, anno_rna = 2, anno_pos = 1, anno_value = 3) {
  check_is_se(ribo)
  check_type(annot, "data.frame", "annot")
  if (!("named_position" %in% colnames(annot))) {
    annot <- .generate_name_positions(annot, anno_rna, anno_pos)
  }

  # Check if annot has the same RNA as the object
  anno_rna_names <- unique(annot[[anno_rna]])
  ribo_rna_names <- S4Vectors::metadata(ribo)$rna_names$current_name
  if (sum((anno_rna_names %in% ribo_rna_names)) == 0) {
    cli::cli_abort(c("Total mismatch in RNA names between annotation and your object !",
      "i" = "Object RNA names : {.val {ribo_rna_names}}.",
      "i" = "Annotation RNA names : {.val {anno_rna_names}}.",
      ">" = "Rename the RNA in the annotation or the object.",
      " " = "(To rename RNA in a object, use {.fn rename_rna})."
    ))
  }

  # Get rowData
  rd <- SummarizedExperiment::rowData(ribo)

  # Construct named_position for checking
  # Assuming rna and rnapos exist
  existing_positions <- paste(rd$rna, formatC(rd$rnapos, width = 4, flag = "0"), sep = "_")

  if (!all(annot[["named_position"]] %in% existing_positions)) {
    missing_positions <- annot[[anno_value]][which(annot[["named_position"]] %in%
      existing_positions == FALSE)]
    len_missing <- length(missing_positions)

    cli::cli_warn(c("{len_missing} position{?s} in your annotation {?is/are} missing in your object !",
      "x" = "Missing positions : {.val {missing_positions}}."
    ))
  }

  # Update site column
  # We match named_position from annot to existing_positions
  # Initialize/Reset site column
  if (!"site" %in% names(rd)) rd$site <- NA

  matches <- match(existing_positions, annot[["named_position"]])
  # If match found, assign anno_value.
  # Note: logic in original: x["site"] <- annot[[anno_value]][match...]

  rd$site <- annot[[anno_value]][matches]

  SummarizedExperiment::rowData(ribo) <- rd

  return(ribo)
}

#' Remove site annotations of a given SummarizedExperiment
#'
#' @param ribo A SummarizedExperiment object.
#' @param annotation_to_remove Specific annotated sites to remove. If set to NULL,
#' the function removes all annotations.
#'
#' @return A SummarizedExperiment object where all annotated position have been replaced by
#' NA, the default value.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' # remove_annotation(ribo_toy)
remove_annotation <- function(ribo, annotation_to_remove = NULL) {
  check_is_se(ribo)
  check_type(annotation_to_remove, "character", "annotation_to_remove")
  rd <- SummarizedExperiment::rowData(ribo)

  if (is.null(annotation_to_remove)) {
    rd$site <- NA
  } else {
    current_annotation <- get_annotation(ribo)
    if (nrow(current_annotation) > 0) {
      new_annotation <- current_annotation[!(current_annotation[["site"]] %in% annotation_to_remove), ]
      # This logic seems circular in original too. It calls keep_selected_annotation.
      # keep_selected calls remove then annotate.
      # simpler:

      # Just set site to NA where site is in annotation_to_remove
      rd$site[rd$site %in% annotation_to_remove] <- NA
    }
  }

  SummarizedExperiment::rowData(ribo) <- rd
  return(ribo)
}

#' Keep only a subset of the current annotation
#'
#' @param ribo a SummarizedExperiment
#' @param annotation_to_keep vector containing annotated sites'name to keep
#'
#' @return a SummarizedExperiment where only annotated sites within annotation_to_keep are
#' still annotated.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' # ribo_toy <- rename_rna(ribo_toy)
#' # ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' # ribo_toy <- keep_selected_annotation(ribo_toy, c("28S_Am1310","28S_Cm2848"))
#'
keep_selected_annotation <- function(ribo, annotation_to_keep) {
  check_is_se(ribo)
  check_type(annotation_to_keep, "character", "annotation_to_keep")
  current_annotation <- get_annotation(ribo)
  ribo <- remove_annotation(ribo)

  # Filter annotation df
  new_annotation <- current_annotation[
    which(current_annotation[["site"]] %in% annotation_to_keep),
  ]

  if (nrow(new_annotation) == 0) {
    cli::cli_abort(c("No currently annotated sites matches with your subset !",
      "i" = "Currently annotated sites : {.val {current_annotation[['site']]}}.",
      "i" = "Your subset : {.val {annotation_to_keep}}."
    ))
  }

  # Use the filtered annotation to annotate again
  # Note: get_annotation returns rna, rnapos, site.
  # We need named_position for annotate_site OR passing rna/pos cols.
  # .generate_name_positions handles rna/pos cols.

  # Ensure column names map correctly.
  # get_annotation returns "rnapos","rna","site"

  new_annotation <- .generate_name_positions(new_annotation, "rna", "rnapos")

  # Calling annotate_site. anno_value=3 which is "site" in get_annotation output?
  # get_annotation returns 3 cols. rnapos(1), rna(2), site(3).
  # So anno_rna=2, anno_pos=1, anno_value=3 works.

  ribo <- annotate_site(ribo, new_annotation, anno_rna = 2, anno_pos = 1, anno_value = 3)

  return(ribo)
}

#' Get annotation of a SummarizedExperiment
#'
#' @param ribo A SummarizedExperiment
#'
#' @return A dataframe with the rna name, the position on the rna and the annotated site.
#' @export
#'
#' @examples
#' data("ribo_toy")
#' data("human_methylated")
#' # ribo_toy <- rename_rna(ribo_toy)
#' # ribo_toy <- annotate_site(ribo_toy,human_methylated)
#' # get_annotation(ribo_toy)
get_annotation <- function(ribo) {
  check_is_se(ribo)
  rd <- SummarizedExperiment::rowData(ribo)

  if (!"site" %in% names(rd)) {
    return(data.frame())
  }

  keep <- !is.na(rd$site)
  current_annotation <- as.data.frame(rd[keep, c("rnapos", "rna", "site")])

  rownames(current_annotation) <- NULL

  return(current_annotation)
}
