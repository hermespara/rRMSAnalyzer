#' Check if a vector of strings correspond to valid metadata in a given SummarizedExperiment
#'
#' This function returns nothing if each string matches with a metadata.
#' It will stop and display a formatted error message if one or more strings in
#' the vector are not valid metadata.
#' @param ribo A SummarizedExperiment.
#' @param metadata_name The vector of string to check against ribo's metadata.
#' @keywords internal
#'
check_metadata <- function(ribo, metadata_name) {
  # Get metadata columns that do not exist in ribo's metadata.
  unmatched_elts <- metadata_name[
    which(!(metadata_name %in% names(SummarizedExperiment::colData(ribo))))
  ]
  len_unmatched <- length(unmatched_elts)
  if (len_unmatched > 0) {
    cli::cli_abort(c(
      "You have supplied names that are not part of the object's metadata (colData).",
      "i" = "The object has the following metadata:
      {.val {names(SummarizedExperiment::colData(ribo))}}.",
      "x" = "{len_unmatched} supplied name{?s} {?is/are} not part of the metadata:
      {.val {unmatched_elts}}."
    ))
  }
}
#' Check if a vector of strings correspond to valid samplenames in a given SummarizedExperiment
#'
#' This function returns nothing if each string matches with the samplenames
#' It will stop and display a formatted error message if one or more strings in
#' the vector are not valid samplenames
#' @param ribo A SummarizedExperiment.
#' @param sample_names The vector of string to check against ribo's samplenames
#' @keywords internal
#'
check_sample <- function(ribo, sample_names) {
  unmatched_elts <- sample_names[
    which(!(sample_names %in% colnames(ribo)))
  ]
  len_unmatched <- length(unmatched_elts)

  if (len_unmatched > 0) {
    cli::cli_abort(c(
      "You have supplied samplenames that are not part of the supplied Object",
      "i" = "Sample names in Object : {.val {colnames(ribo)}}.",
      "x" = "{len_unmatched} supplied name{?s} {?is/are} not part of the Object's samplenames:
      {.val {unmatched_elts}}."
    ))
  }
}

#' Check if an object is a SummarizedExperiment
#'
#' @param object The object to check.
#' @keywords internal
check_is_se <- function(object) {
  if (missing(object)) {
    cli::cli_abort("Argument {.arg ribo} is missing, with no default.")
  }
  if (!inherits(object, "SummarizedExperiment")) {
    cli::cli_abort(c(
      "{.arg ribo} must be a SummarizedExperiment object.",
      "x" = "You supplied a {.cls {class(object)}}."
    ))
  }
}

#' Check the type of an argument
#'
#' @param x The argument to check.
#' @param type Expected type ("character", "numeric", "logical", "data.frame", "list").
#' @param name Name of the argument (for error message).
#' @param length Expected length (optional).
#' @keywords internal
check_type <- function(x, type, name, length = NULL) {
  if (missing(x) || is.null(x)) {
    return(invisible(NULL))
  } # Handle optional args separately if needed, or assume caller handles NULL

  valid <- switch(type,
    "character" = is.character(x),
    "numeric" = is.numeric(x),
    "logical" = is.logical(x),
    "data.frame" = is.data.frame(x),
    "list" = is.list(x),
    "matrix" = is.matrix(x),
    FALSE
  )

  if (!valid) {
    cli::cli_abort("{.arg {name}} must be of type {.cls {type}}.")
  }

  if (!is.null(length)) {
    if (length(x) != length) {
      cli::cli_abort("{.arg {name}} must have length {length}.")
    }
  }
}

#' Check if an argument is in a allowed set of values
#'
#' @param x The argument to check.
#' @param set The set of allowed values.
#' @param name Name of the argument (for error message).
#' @keywords internal
check_in_set <- function(x, set, name) {
  if (missing(x) || is.null(x)) {
    return(invisible(NULL))
  }

  if (!all(x %in% set)) {
    cli::cli_abort(c(
      "{.arg {name}} must be one of {.val {set}}.",
      "x" = "You supplied {.val {x}}."
    ))
  }
}

#' Validate a named color mapping against observed metadata values.
#'
#' @param colors Named character vector of colors.
#' @param values Metadata values that will be mapped to colors.
#' @param arg_name Argument name used in error messages.
#' @keywords internal
check_named_colors <- function(colors, values, arg_name = "sample_colors") {
  if (missing(colors) || is.null(colors)) {
    return(invisible(NULL))
  }

  check_type(colors, "character", arg_name)

  if (is.null(names(colors)) || any(names(colors) == "")) {
    cli::cli_abort("{.arg {arg_name}} must be a named character vector.")
  }

  values <- unique(as.character(values[!is.na(values)]))
  missing_values <- setdiff(values, names(colors))

  if (length(missing_values) > 0) {
    cli::cli_abort(c(
      "{.arg {arg_name}} must provide a color for each plotted group.",
      "x" = "Missing color mapping for {.val {missing_values}}."
    ))
  }

  invisible(NULL)
}

# Compute ellipse coordinates and centroids for grouped 2D data.
.group_ellipses <- function(data, x_col, y_col, group_col,
                            level = 0.95, npoints = 100,
                            regularization = 1e-4) {
  if (is.null(group_col) || !group_col %in% names(data)) {
    return(list(ellipses = NULL, centroids = NULL))
  }

  plot_df <- data[, c(x_col, y_col, group_col), drop = FALSE]
  colnames(plot_df) <- c("x", "y", "group")
  plot_df <- stats::na.omit(plot_df)
  plot_df$group <- as.character(plot_df$group)

  if (!nrow(plot_df)) {
    return(list(ellipses = NULL, centroids = NULL))
  }

  centroids <- stats::aggregate(cbind(x, y) ~ group, data = plot_df, FUN = mean)
  group_split <- split(plot_df[, c("x", "y")], plot_df$group)

  ellipses <- lapply(names(group_split), function(group_name) {
    group_df <- group_split[[group_name]]
    if (nrow(group_df) < 3) {
      return(NULL)
    }

    cov_matrix <- stats::cov(group_df)
    if (any(!is.finite(cov_matrix))) {
      return(NULL)
    }

    diag_scale <- max(mean(diag(cov_matrix)), .Machine$double.eps)
    cov_matrix <- cov_matrix + diag(regularization * diag_scale, 2)
    eig <- eigen(cov_matrix, symmetric = TRUE)
    eig$values[eig$values < .Machine$double.eps] <- .Machine$double.eps

    angles <- seq(0, 2 * pi, length.out = npoints)
    unit_circle <- rbind(cos(angles), sin(angles))
    radius <- sqrt(stats::qchisq(level, df = 2))
    shape <- eig$vectors %*% diag(sqrt(eig$values), nrow = 2) %*% unit_circle * radius
    center <- c(mean(group_df$x), mean(group_df$y))

    data.frame(
      x = center[1] + shape[1, ],
      y = center[2] + shape[2, ],
      group = group_name
    )
  })

  ellipses <- dplyr::bind_rows(ellipses)
  if (!nrow(ellipses)) {
    ellipses <- NULL
  }

  return(list(ellipses = ellipses, centroids = centroids))
}
