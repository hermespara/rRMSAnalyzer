#' Create a RiboClass from count files and metadata.
#'
#' @param count_path (required) path to the data folder containing count files.
#' @param metadata  Data frame or path to a CSV file containing metadata.
#' @param count_sep Delimiter used for the CSV files.
#' @param metadata_sep Delimiter used in metadata (if imported from file).
#' @param count_header Boolean, specify if count files have a header or not.
#' @param count_value Name or index of the column containing count values.
#' @param count_rnaid Name or index of the column containing the name of the RNA in count data.
#' @param count_pos Name or index of the column containing the site's position in count data.
#' @param metadata_key (required when metadata not null) Name or index of the column containing the samples' filename.
#' @param metadata_id Name or index of the column containing the sample name.
#'
#' @description
#' Read ribomethseq count files and their associated metadata and merge them into a SummarizedExperiment.
#'
#' @details
#' The resulting \code{SummarizedExperiment} object contains:
#'
#' * \code{assays}: A list containing the \code{counts} matrix.
#' * \code{rowData}: A \code{DataFrame} containing information for each site:
#'   - \code{rna}: The name of the RNA.
#'   - \code{rnapos}: The position on the RNA.
#'   - \code{site}: The site annotation (if applicable).
#' * \code{colData}: A \code{DataFrame} containing metadata for each sample.
#' * \code{metadata}: A list containing \code{rna_names} table.
#'
#' The path given in \code{count_path} should contain only necessary count files (one per sample).
#' While the directory structure is not important, make sure each sample has a unique filename.
#'
#' The \code{metadata} must contain a column (specified by \code{metadata_key}) that matches the
#' filenames of the count files.
#' @md
#'
#' @return A SummarizedExperiment object.
#'
create_se <- function(count_path,
                      metadata = NULL,
                      count_sep = "\t",
                      metadata_sep = ",",
                      count_header = FALSE,
                      count_value = 3,
                      count_rnaid = 1,
                      count_pos = 2,
                      metadata_key = "filename",
                      metadata_id = NULL) {
  # Parameter validation
  check_type(count_path, "character", "count_path", length = 1)
  check_type(count_sep, "character", "count_sep", length = 1)
  check_type(metadata_sep, "character", "metadata_sep", length = 1)
  check_type(count_header, "logical", "count_header", length = 1)
  check_type(metadata_key, "character", "metadata_key", length = 1)

  if (!is.null(metadata_id) && !is.character(metadata_id) && !is.numeric(metadata_id)) {
    cli::cli_abort("{.arg metadata_id} must be a character string or a numeric index.")
  }

  # loading metadata
  if (is.null(metadata)) {
    metadata <- generate_metadata_df(count_path, create_samplename_col = FALSE)
    rna_counts_dt <- .read_count_files(count_path, count_sep, count_header, count_rnaid, count_pos, count_value)
    if (length(rna_counts_dt) == 0) cli::cli_abort("ERROR : no file was loaded")
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])
  } else {
    if (is.character(metadata)) {
      if (file.exists(metadata)) {
        metadata <- utils::read.csv(metadata, sep = metadata_sep)
      } else {
        cli::cli_abort("the path specified for metadata does not exist or is not a file !")
      }
    }

    if (!is.data.frame(metadata)) {
      cli::cli_abort("metadata must be a dataframe or a path to a csv file !")
    }

    if (is.character(metadata_key) && !(metadata_key %in% names(metadata))) {
      cli::cli_abort("{metadata_key} (metadata_key param) is not a column in metadata")
    }

    if (is.character(metadata_id) && !(metadata_id %in% names(metadata))) {
      cli::cli_abort("{metadata_id} (metadata_id param) is not a column in metadata")
    }

    # rename the column specified in "metadata_id" to "samplename"
    if (is.character(metadata_id)) {
      names(metadata)[names(metadata) == metadata_id] <- "samplename"
    } else {
      names(metadata)[metadata_id] <- "samplename"
    }

    # Check if metadata has not duplicated filename or samplename
    if (anyDuplicated(metadata[metadata_key])) {
      cli::cli_abort("ERROR! Duplicated filename in metadata. Each sample must have an unique filename.")
    }

    if (anyDuplicated(metadata["samplename"])) {
      cli::cli_abort("ERROR! Duplicated samplename in metadata. Each sample must have an unique samplename.")
    }

    rownames(metadata) <- metadata[, "samplename"]

    # read count data
    rna_counts_dt <- .read_count_files(count_path, count_sep, count_header,
      count_rnaid, count_pos, count_value,
      metadata_filenames = as.character(metadata[, metadata_key])
    )


    # generate RNA names table
    rna_names_df <- .generate_rna_names_table(rna_counts_dt[[1]])

    # keep only metadatas that have associated data

    metadata <- metadata[which(metadata[, metadata_key] %in% names(rna_counts_dt)), ]

    # Rename sample in counts list according to the names in metadata
    names(rna_counts_dt) <- metadata[, "samplename"][match(names(rna_counts_dt), metadata[, metadata_key])]

    # order samples by metadata and remove null samples
    rna_counts_dt <- rna_counts_dt[metadata[, "samplename"]]
    rna_counts_dt <- rna_counts_dt[lengths(rna_counts_dt) != 0]
  }

  # Construction of SummarizedExperiment

  # 1. Assay (Counts)
  # Extract counts from each dataframe to form a matrix
  # We assume all dataframes have the same positions (checked in .read_count_files)
  counts_matrix <- vapply(rna_counts_dt, function(x) x$count, numeric(nrow(rna_counts_dt[[1]])))

  # 2. RowData (Positions)
  # Use the first sample to define positions
  row_data <- rna_counts_dt[[1]][, c("rna", "rnapos", "site")]
  # Make sure 'rna' is a factor as expected downstream
  row_data$rna <- factor(row_data$rna)

  # 3. ColData (Metadata)
  col_data <- metadata

  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowData = row_data,
    colData = col_data,
    metadata = list(rna_names = rna_names_df) # Store rna_names in generic metadata slot
  )

  return(se)
}


#' Import and transform a list of count files for a SummarizedExperiment
#' @description
#' This internal function is used to read all count files and transform them into a list of dataframes.
#'
#' @param path_to_files path to the folder containing count files
#' @param sep delimiter used in count files
#' @param rna_col Name or index of the column containing RNA names
#' @param position_col Name or index of the column containing the site positions
#' @param count_value Name or index of the column containing the count values
#' @param header Boolean, specify if count files have a header
#' @param metadata_filenames Filenames in metadata to check against (optional)
#'
#' @return a list of sample dataframes
#'
#' @keywords internal
.read_count_files <-
  function(path_to_files,
           sep,
           header,
           rna_col,
           position_col,
           count_value,
           metadata_filenames = NULL) {
    if (!dir.exists(path_to_files)) cli::cli_abort("the path given for the csv files does not exist or is not a directory !")

    rna_counts_fl <-
      list.files(path_to_files, recursive = TRUE, full.names = TRUE)

    # Check if the filenames on disk match filenames in metadata. Fail otherwise.
    if (!is.null(metadata_filenames)) {
      pat <- paste0("\\b(", paste(metadata_filenames, collapse = "|"), ")\\b")
      rna_counts_fl <- rna_counts_fl[grep(pat, rna_counts_fl)]

      if (length(metadata_filenames) > length(rna_counts_fl)) {
        rna_count_missing <- metadata_filenames[which(
          !is.element(metadata_filenames, basename(rna_counts_fl))
        )]
        cli::cli_warn(paste("File", rna_count_missing, "does not exist. Typo in metadata ?", collapse = "\n"))
      }
    }

    # 1) Check if there is any duplicated name in the filename list

    if (anyDuplicated(basename(rna_counts_fl)) > 0) {
      cli::cli_abort("ERROR: some samples share the same filename!")
    }

    # Check if there is a mismatch between the filenames in metadata
    # and the filenames in the folder

    rna_counts_dt <-
      lapply(rna_counts_fl, utils::read.csv, sep = sep, header = header)

    # stop if total mismatch between filenames and the filename column given in metadata
    if (length(rna_counts_dt) == 0) cli::cli_abort(paste0("ERROR! No file has the filenames specified in the metadata column specified by metadata_key"))

    # check if there are less than 3 columns, which can happen when one fails to
    # specify the correct separator
    first_file <- rna_counts_dt[[1]]
    if (ncol(first_file) < 3) {
      cli::cli_abort("not enough columns in your count data !
      \nCheck if you have specified the correct columns separator in count_sep")
    }

    # Helper function to check column existence
    check_col <- function(col, df, arg_name) {
      if (is.character(col)) {
        if (!(col %in% names(df))) {
          cli::cli_abort(c(
            "x" = "{.arg {col}} ({arg_name}) is not present in data",
            "i" = "Please check your SummarizedExperiment initialization parameters",
            "*" = "Valid columns are: {.val {colnames(df)}}"
          ))
        }
      } else if (is.numeric(col)) {
        if (col > ncol(df) || col < 1) {
          cli::cli_abort("Column index {.val {col}} ({arg_name}) is out of bounds (file has {ncol(df)} columns).")
        }
      } else {
        cli::cli_abort("{arg_name} must be a character string or numeric index.")
      }
    }

    check_col(rna_col, first_file, "count_rnaid")
    check_col(position_col, first_file, "count_pos")
    check_col(count_value, first_file, "count_value")

    # 2) reorder cols for each count table and name them like this :
    # RNA | position_on_rna | count
    # and add a siteID column with a default value of NA
    rna_counts_dt <- lapply(rna_counts_dt, function(x) {
      x <- x[, c(rna_col, position_col, count_value)]
      colnames(x) <- c("rna", "rnapos", "count")
      x["site"] <- NA
      return(x)
    })

    # 3) combine both the  # re-calculate named_position column for each samples.
    # This column name 'named_position' is hardcoded in extract_data()
    rna_counts_dt <- .generate_se_named_position(rna_counts_dt, "rna", "rnapos")


    # 4) give a name for each element of the list
    # TODO: Check if the elements in RNA_counts_dt are in the same order as in RNA_counts_fl
    names(rna_counts_dt) <- basename(rna_counts_fl)

    # 5) using named_position, we check if samples share the same positions
    reference_sample_name <- names(rna_counts_dt)[1]

    sample_check_results <-
      vapply(names(rna_counts_dt), function(x) {
        .check_sample_positions(
          sample_1 = rna_counts_dt[[1]],
          sample_1_name = reference_sample_name,
          sample_2 = rna_counts_dt[[x]],
          sample_2_name = x
        )
      }, logical(1))

    failing_samples <-
      names(sample_check_results[which(sample_check_results == FALSE)])
    if (identical(failing_samples, character(0))) {
    } else {
      cli::cli_abort(paste(
        "[ERROR] The following samples have failed the positions check : ",
        paste(failing_samples, collapse = "; ")
      ))
    }



    return(rna_counts_dt)
  }
