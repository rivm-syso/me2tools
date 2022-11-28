#' Get the values contributing to the Qmain from the results file (.txt)
#'
#' Typical use of Multilinear Engine version 2 (ME-2) will provide three
#' different files after a run: \dQuote{.dat}: machine readable format,
#' \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file with
#' auxiliary information (i.e. headers). The headers in the text (.txt) file
#' are used in me2tools to split the data into several blocks. In this function
#' the values contributing to Qmain are read using user provided header
#' information, denoting the \dQuote{start} and \dQuote{end} lines of the block
#' containing the required data.
#'
#' @param me2_txt_file ME2 output file (.txt), containing the results and
#'   auxiliary information.
#'
#' @return tibble containing Qmain and contributions to Qmain for multiple runs
#'
#' @export
#'
#' @seealso \code{\link{me2_read_all}}
#'
#' @import cli
#' @import readr
#' @import stringr
#' @import tibble
#' @import purrr
#' @import dplyr
#'
me2_read_Qcontrib <- function(me2_txt_file,
                              block_boundaries = list(
                                "start" = "^\\s*Contribution in Qmain by different variables",
                                "end" = "^\\s*Task #,  seed:"
                              )) {

  # check if file exists
  if (!file.exists(me2_txt_file)) {
    cli::cli_abort(c(
      "{.var me2_txt_file} must be a valid file:",
      "x" = "File '{me2_txt_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_txt_file)
  text_filter <- text[text != ""]

  # check if we only have one parameter block
  me2_dupl_block_check(text_filter)

  # factor profiles are located between "   Contribution in Qmain by different
  # variables" and "    Task #,  seed:    x  x" (note the spaces, and with x
  # as numbers)
  index_start <- stringr::str_which(text_filter, block_boundaries$start)
  if (length(index_start) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The start header of the data block could not be found in '{file}'.",
      "x" = "Did you provide the correct file?"
    ))
  }

  index_end <- stringr::str_which(text_filter, block_boundaries$end)

  # index_end now contains the first line and not the last one.
  index_end <- c(index_end[2:length(index_end)], length(text_filter))
  if (length(index_end) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The end header of the data block could not be found in '{file}'.",
      "x" = "Did you provide the correct file?"
    ))
  }

  if (length(index_start) != length(index_end)) {
    if(length(index_end) == 2*length(index_start)) {
      # replace index_end with the even rows
      row_odd <- seq_len(length(index_end)) %% 2
      index_end <- index_end[row_odd == 0]
    } else {
      cli::cli_abort(c(
        "Length of start and end headers do not match:",
        "i" = "The lengths of found start and end headers are not equal",
        "x" = "Does the file contain duplicate entries within a block of data?"
      ))
    }
  }

  # the length of the index_start will provide how many runs were performed.
  # we need to get all runs
  contrib_matrix <- tibble::tibble()
  for (run.number in seq(1, length(index_start), 1)) {
    matrix.text <- text_filter[seq(index_start[run.number] + 1, index_end[run.number] - 1, 1)]
    contrib_matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
      delim = " ",
      trim_ws = TRUE,
      col_names = FALSE
    )) %>%
      purrr::set_names(c("species", "contrib.Qmain")) %>%
      tibble::add_column(model_type = "ME-2", .before = "species") %>%
      tibble::add_column(model_run = run.number, .before = "species")


    if (nrow(contrib_matrix) > 0) {
      contrib_matrix <- dplyr::bind_rows(contrib_matrix, contrib_matrix.tmp)
    } else {
      contrib_matrix <- contrib_matrix.tmp
    }
    rm(contrib_matrix.tmp)
  }
  browser()
  return(contrib_matrix)
}
