#' Extract ME2 matrix between two line numbers
#'
#' This function is used to extract text from between two line numbers
#' previously determined. This is an internal function.
#'
#' @param text text
#' @param line.start the position in the text where the line starts
#' @param line.end the position in the text where the line ends
#' @param headers when set to \code{TRUE} the default column headers are
#'   generated based on the available information in the text. If set to
#'   \code{FALSE} this function will provide default headers.
#'
#' @noRd
#'
#' @import readxl
#' @import dplyr
#' @import lubridate
#' @import stringr
#'
extract_me2_matrix_between <- function(text,
                                       line.start,
                                       line.end,
                                       headers = TRUE) {
  text_filter <- text[text != ""]

  if (headers) {
    # grab column names and make them compatible to PMFR
    colnames <- strsplit(text_filter[line.start + 1], "\\s+")
    colnames[[1]] <- tolower(colnames[[1]])
    colnames[[1]][1] <- "identifier"
    # set factor names similar to EPA-PMF
    # colnames[[1]] <- gsub("^0", "" , colnames[[1]])
    colnames[[1]] <- gsub("fact", "factor_", colnames[[1]])
    matrix.text <- text_filter[seq(line.start + 2, line.end - 1, 1)]

    # concentration of species
    matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
      delim = " ",
      trim_ws = TRUE,
      col_names = FALSE
    )) %>%
      select(-dplyr::first(names(.)), -dplyr::last(names(.))) %>%
      purrr::set_names(colnames[[1]])
  } else {
    matrix.text <- text_filter[seq(line.start + 1, line.end - 1, 1)]

    # concentration of species
    matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
      delim = " ",
      trim_ws = TRUE,
      col_names = FALSE
    )) %>%
      select(-dplyr::first(names(.))) # drop the first column with only numbers

    colnames <- paste0("factor_", sprintf("%02d", seq(1, ncol(matrix.tmp), 1)))

    matrix.tmp <- matrix.tmp %>%
      purrr::set_names(colnames)
  }

  return(matrix.tmp)
}
