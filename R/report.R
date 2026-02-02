#' Function to generate HTML report using Rmarkdown
#' @description
#' Note : if you want to create your own report template,
#' please make sure it has the following header :
#'
#' ---
#' title: "RMSAnalyzer : title of the report"
#' output: html_document
#'
#' params:
#'  library_col: "library"
#'  condition_col: "condition"
#'  ribo_name: "ribo"
#'  project_name: "Unnamed project"
#' ---
#'
#' @param ribo a SummarizedExperiment object.
#' @param rmdfile filename of the Rmarkdown template to use.
#' @param condition_col name or index of the column __in metadata__ containing the condition
#' @param library_col name or index of the column __in metadata__ containing the library
#' @param project_name Name of the project to display on the report.
#' @param output_dir Path of output directory
#'
#' @keywords internal
#'
report <- function(ribo, ribo_name, rmdfile, condition_col, library_col, project_name, output_dir) {
  ribo_to_check <- ribo # used in Rmarkdown template # nolint
  path <- system.file("rmd", package = "rRMSAnalyzer")
  rmarkdown::render(file.path(path, rmdfile),
    params = list(
      library_col = library_col,
      condition_col = condition_col,
      ribo_name = ribo_name,
      project_name = project_name
    ),
    output_dir = output_dir
  )
}

#' Generate html report about possible batch effect for a given SummarizedExperiment.
#'
#' @param ribo A SummarizedExperiment object.
#' @param library_col Library/run column in the metadata.
#' @param project_name Name of the project.
#' @param output_dir Path to output dir. Working directory is selected by default.
#' @export
report_qc <- function(ribo, library_col, project_name = "Unnamed project", output_dir = getwd()) {
  check_is_se(ribo)
  check_type(library_col, "character", "library_col", length = 1)
  check_type(project_name, "character", "project_name", length = 1)
  check_type(output_dir, "character", "output_dir", length = 1)
  check_metadata(ribo, library_col)
  ribo_name <- deparse(substitute(ribo))

  report(ribo, ribo_name, "quality_control.Rmd",
    condition_col = NULL,
    library_col = library_col,
    project_name = project_name,
    output_dir = output_dir
  )
}
