#' Get the factor "correlations" from the BS results file
#'
#' Typical use of Bootstrap using Multilinear Engine version 2 (ME-2) will
#' provide three different files after a run: \dQuote{.dat}: machine readable
#' format, \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file
#' with auxiliary information (i.e. headers). The headers in the text (.txt)
#' file are used in me2tools to split the data into several blocks. In this
#' function the factor "correlations" with the Best-fit factors are read using
#' the provided header information, denoting the \dQuote{start} line of the
#' block containing the required data.
#'
#' @param me2_bs_txt_file file name of the BS results file.
#' @param tidy_output Should the output be reshaped into tidy data? Default:
#'   FALSE
#' @param block_boundaries A list containing the \dQuote{start} string used to
#'   identify the boundaries of the block containing the values for the
#'   correlations. Typically, this block is at the end of the file and therefor
#'   no \dQuote{end} string is required.
#'
#' @return tibble containing the factor "correlations" with Best-fit factors
#'
#' @export
#'
#' @import cli
#' @import readr
#' @import stringr
#' @import dplyr
#'
me2_BS_read_correlations <- function(me2_bs_txt_file,
                                     tidy_output = FALSE,
                                     block_boundaries = list("start" = "factor \"correlations\" with Best-fit factors")) {

  # check if file exists
  if (!file.exists(me2_bs_txt_file)) {
    cli::cli_abort(c(
      "{.var me2_bs_txt_file} must be a valid file:",
      "x" = "File '{me2_bs_txt_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_bs_txt_file)
  text_filter <- text[text != ""]

  # factor profiles are located between "factor \"correlations\" with Best-fit
  # factors".
  index_start <- stringr::str_which(text_filter, block_boundaries$start)
  if (length(index_start) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The start header of the data block could not be found in '{me2_bs_txt_file}'.",
      "x" = "Did you provide the correct file or the correct start string using
      {var. block_boundaries}?"
    ))
  }
  # The correlation matrix is the latest block at the end of the file.
  index_end <- length(text_filter)

  matrix.text <- text_filter[seq(index_start + 1, index_end, 1)]

  # concentration of species
  matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
    delim = " ",
    trim_ws = TRUE,
    col_names = FALSE
  ))
  col_names <- c(
    c("model_run", "Qmain", "Qaux", "Q_XX"),
    paste0("factor_", sprintf("%02d", seq(1, ncol(matrix.tmp) - 4, 1)))
  )

  names(matrix.tmp) <- col_names

  if (tidy_output) {
    matrix.tmp <- matrix.tmp %>%
      pivot_longer(cols = !c("model_run", "Qmain", "Qaux", "Q_XX"),
                   names_to = "factor",
                   values_to = "corr")
  }
  return(matrix.tmp)
}
