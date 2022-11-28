#' Check for duplicate outputs in the ME-2 output file
#'
#' Duplicate blocks of data can be outputed at several locations in the ME-2
#' ini file. In order to read the data, only 1 block of data can be present in
#' the output.
#'
#' Code execution is stopped when this is the case.
#' 
#' @param text
#'
#' @noRd
#'
#' @import stringr
#' @import cli

me2_dupl_block_check <- function(text) {
  n_block <- sum(stringr::str_count(text, "Contributions to Qaux by smooth & norm equ:s, for each factor"))
  n_time_factors <- sum(stringr::str_count(text, "Time factors"))

  if (n_block > n_time_factors) {
    cli::cli_abort(c(
      "Can't have duplicate data block:",
      "i" = "The ME-2 output file contains a duplicate block with parameter
        settings'.",
      "x" = "Check your ME-2 ini file (especially the postproc section) for
        duplicate output."
    ))
  }
}
