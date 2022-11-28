#' Get data from the ME-2 log file
#'
#' The Multilinear Engine, version 2 (ME-2), has several options for default
#' logging. This function tries to read the messages from the log file. Please
#' note that inside the .ini file you can set the amount of output for the log
#' file. In some cases these outputs are very verbose and cannot be used in this
#' function. In this function the messages inside the LOG file are read using the
#' following header information, denoting the \dQuote{start} and \dQuote{end}
#' lines of the block containing the required data:
#'
#' @param me2_log_file file name of the ME-2 log file.
#' @param block_boundaries A list containing the \dQuote{start} string used
#'   to identify the boundaries of the block containing the values
#'   for the log file.
#'
#' @return tibble containing the following information over each base_run:
#'   \dQuote{convergence} (T/F); \dQuote{self.cancel} (T/F);
#'   \dQuote{iter.steps} (#); \dQuote{precon_mode} (#) and \dQuote{messages}.
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
me2_read_log <- function(me2_log_file,
                         block_boundaries = list("start" = "^\\s*Starting to process with the same .INI file next task:")) {

  # check if file exists
  if (!file.exists(me2_log_file)) {
    cli::cli_abort(c(
      "{.var me2_log_file} must be a valid file:",
      "x" = "File '{me2_log_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_log_file)
  text_filter <- text[text != ""]

  # log message are located between "Starting to process with the same .INI file next task:"
  index_start <- stringr::str_which(text_filter, block_boundaries$start)
  if (length(index_start) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The start header of the data block could not be found in '{file}'.",
      "x" = "Did you provide the correct file, or did you limit the output
        written to the log file in the .ini file?"
    ))
  }
  index_start <- c(1, index_start)
  index_end <- c(index_start[2:length(index_start)], length(text_filter))
  log_me2 <- tibble::tibble()
  for (run.number in seq(1, length(index_start), 1)) {
    run.text <- text_filter[seq(index_start[run.number], index_end[run.number] - 1, 1)]

    # search for no-convergence message
    convergence <- TRUE
    me2.message <- NULL
    index <- stringr::str_which(run.text, "Iteration interrupted because maximum step count")
    if (length(index) == 1) {
      me2.message <- gsub("\\*\\\\  ", "", gsub("\r?\n|\r?^\\s*\\*/\\s*", "", paste0(run.text[index], run.text[index + 1])))
      convergence <- FALSE
    }

    # search for self-cancel limit
    self_cancel <- FALSE
    index <- stringr::str_which(run.text, "Iteration interrupted, exceeding self-cancel limit")
    if (length(index) == 1) {
      tmp.me2.message <- gsub("\\*\\\\  ", "", gsub("\r?\n|\r?^\\s*\\*/\\s*", "", paste0(run.text[index], run.text[index + 1])))
      if (identical(me2.message, NULL)) {
        me2.message <- tmp.me2.message
      } else {
        me2.message <- paste(me2.message,
          tmp.me2.message,
          sep = ";"
        )
      }
      self_cancel <- TRUE
    }
    # replace null with NA
    if(is.null(me2.message)) {
      me2.message <- NA
    }

    # get iter-steps
    line <- "iter-steps in all in task number\\s*[0-9]*$"
    if (length(stringr::str_which(run.text, line)) > 0) {
      iter_steps <- as.numeric(gsub(line, "", run.text[stringr::str_which(run.text, line)]))
    } else {
      iter_steps <- NA
    }


    # get iter-steps
    line <- "Preconditioning mode =\\s*"
    if (length(stringr::str_which(run.text, line)) > 0) {
      precon_mode <- as.numeric(gsub(line, "", run.text[stringr::str_which(run.text, line)]))
    } else {
      precon_mode <- NA
    }

    # create tibble
    log_tmp <- tibble::tibble(
      "base_run" = run.number,
      "convergence" = convergence,
      "self.cancel" = self_cancel,
      "iter.steps" = iter_steps,
      "precon.mode" = precon_mode,
      "message" = me2.message
    )

    if (nrow(log_me2) > 0) {
      log_me2 <- dplyr::bind_rows(log_me2, log_tmp)
    } else {
      log_me2 <- log_tmp
    }
    rm(log_tmp)
  }
  return(log_me2)
}
