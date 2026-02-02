#' Known 2'Ome positions in Humans' rRNA
#'
#' A dataset containing 112 2'Ome positions with associated data
#'
#' @usage data(human_methylated)
#'
#' @format A data frame with 112 rows and 10 variables:
#' \describe{
#'   \item{Position}{nucleotide position on rRNA}
#'   \item{rRNA}{rRNA where this methylation is found}
#'   \item{Nomenclature}{name given to this modification}
#'   \item{NR_046235.Numbering}{Position on}
#'   \item{Sequence}{nucleotides sequence surround position}
#'   \item{SNORD}{ID of associated snoRNA}
#'   \item{Mode.of.coding}{how the SNORD is coded}
#'   \item{SNORD.host.gene}{SNORD's host gene}
#'   \item{Ensembl}{SNORD's Ensembl reference}
#'   \item{Nucleotide}{the nucleotide present at the position (A,T,G or C)}
#'   ...
#' }
"human_methylated"

#' Suspected 2'Ome positions in Humans' rRNA
#'
#' A dataset containing 17 2'Ome positions with associated data
#' @usage data(human_suspected)
#'
#' @format A data frame with 17 rows and 10 variables:
#' \describe{
#'   \item{Position}{nucleotide position on rRNA}
#'   \item{rRNA}{rRNA where this methylation is found}
#'   \item{Nomenclature}{name given to this modification}
#'   \item{NR_046235.Numbering}{Position on}
#'   \item{Sequence}{nucleotides sequence surround position}
#'   \item{SNORD}{ID of associated snoRNA}
#'   \item{Mode.of.coding}{how the SNORD is coded}
#'   \item{SNORD.host.gene}{SNORD's host gene}
#'   \item{Ensembl}{SNORD's Ensembl reference}
#'   \item{Nucleotide}{the nucleotide present at the position (A,T,G or C)}
#'   ...
#' }
"human_suspected"

#' SummarizedExperiment from a toy dataset
#'
#' A SummarizedExperiment object containing 10 samples + 2 reference RNA.
#'
#' Samples are from 4 different biological conditions ("condition" column in metadata).
#' The sequencing has been done in two different batches ("run" column in metadata).
#' Both batches have the same reference RNA, to detect technical bias.
#' @usage data(ribo_toy)
#'
#' @format a SummarizedExperiment object with:
#' \describe{
#'   \item{assays}{list containing counts and cscore matrices}
#'   \item{colData}{metadata dataframe of all samples}
#'   \item{rowData}{dataframe containing positions, RNA names and site annotations}
#'   \item{metadata}{list containing the rna_names table}
#'   }
"ribo_toy"
