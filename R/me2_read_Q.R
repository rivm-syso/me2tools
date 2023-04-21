#' Get the Q values from the results file (.txt)
#'
#' Typical use of Multilinear Engine version 2 (ME-2) will provide three
#' different files after a run: \dQuote{.dat}: machine readable format,
#' \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file with
#' auxiliary information (i.e. headers). The headers in the text (.txt) file
#' are used in me2tools to split the data into several blocks. In this function
#' the Q values are read using user provided header information, denoting the
#' \dQuote{start} line of the block containing the required data.
#' 
#' This function can also read the Q values created by the modified EPA-PMF
#' ME2 script. While reading the file, the function checks for both formats
#' prior to using the \code{block_boundaries} and tries to read the Q values 
#' based on the format found. The \code{block_boundaries} are not used when 
#' reading data from the modified EPA-PMF ME2 script.
#'
#' @param me2_txt_file ME2 output file (.txt), containing the results and
#'   auxiliary information.
#' @param block_boundaries A list containing the \dQuote{start} string used to
#'   identify the boundaries of the block containing the Q values. The 
#'   \dQuote{end} string is automatically derived as the Q values is typically
#'   the last block in the file.
#'
#' @return tibble containing the Q values for each run.
#'
#' @export
#'
#' @seealso \code{\link{me2_read_all}}
#'
#' @import cli
#' @import readr
#' @import stringr
#' @import tibble
#' @import dplyr
#'
me2_read_Q <- function(me2_txt_file,
                       block_boundaries = list("start" = "^\\s*Sum-of-squares  Q, Qmain, Qaux")) {

  # check if file exists
  if (!file.exists(me2_txt_file)) {
    cli::cli_abort(c(
      "{.var me2_txt_file} must be a valid file:",
      "x" = "File '{me2_txt_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_txt_file)
  text_filter <- text[text != ""]

  # The file can also be from the modified EPA script, so if this is the case
  # we need to read the Q-values in an other way.
  
  if(stringr::str_detect(text_filter[[1]], "iniparams.txt")) {
    # this is likely from the EPA modified script.
    q_values <- me2_BS_read_correlations(me2_bs_txt_file = me2_txt_file) %>% 
      dplyr::select(- tidyr::all_of(tidyr::contains("factor")))
    
    names(q_values)[1] <- "model_run"
    
    cli::cli_alert_info(c(
      "Input looks like modified EPA-script, returning found Q-values."
    ))
    
  } else {
    # check if we only have one parameter block
    me2_dupl_block_check(text_filter)
  
    # Get the Q-values (at the line containing "  Sum-of-squares  Q, Qmain, Qaux")
    index_start <- stringr::str_which(text_filter, block_boundaries$start)
  
    if (length(index_start) == 0) {
      cli::cli_abort(c(
        "Header not found:",
        "i" = "The start header of the data block could not be found in '{me2_txt_file}'.",
        "x" = "Did you provide the correct file?"
      ))
    }
  
    q_values <- tibble::tibble()
    for (run.number in seq(1, length(index_start), 1)) {
      split_result <- strsplit(
        trimws(text_filter[index_start[run.number]],
          which = "both"
        ),
        "\\s+"
      )
  
      q_values_tmp <- tibble::tibble(
        model_run = run.number,
        Q = as.numeric(split_result[[1]][5]),
        Qmain = as.numeric(split_result[[1]][6]),
        Qaux = as.numeric(split_result[[1]][7])
      )
  
      if (nrow(q_values) > 0) {
        q_values <- dplyr::bind_rows(q_values, q_values_tmp)
      } else {
        q_values <- q_values_tmp
      }
      rm(q_values_tmp)
    }
  }

  return(q_values)
}
