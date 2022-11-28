#' Get the residuals from the residual file
#'
#' Typical use of Multilinear Engine version 2 (ME-2) will provide three
#' different files after a run: \dQuote{.dat}: machine readable format,
#' \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file with
#' auxiliary information (i.e. headers). The headers in file containing the
#' residual results (.rsd) file are used in me2tools to split the data into
#' several blocks. In this function the residuals are read using user provided
#' header information, denoting the \dQuote{start}line of the block containing
#' the required data.
#'
#' @param me2_res_file ME2 output file (.rsd), containing the results for the
#'   residuals.
#' @param block_boundaries A list containing the \dQuote{start} string used to
#'   identify the boundaries of the block containing the values for the
#'   residuals. For residuals the \dQuote{end} string is automatically derived
#'   based on the \dQuote{start} string for the next block.
#' @param species A vector containing the names of the species for the columns
#'   in the residuals.
#' @param dates A vector containing the sample dates for the rows in the
#'   residuals.
#' @param tz Parameter to control the time zone when parameter \dQuote{dates} is
#'   not used. Default: 'Etc/GMT-1'
#'
#' @return tibble containing scaled residuals for each run.
#'
#' @export
#'
#' @seealso \code{\link{me2_read_all}}
#'
#' @import cli
#' @import readr
#' @import stringr
#' @import dplyr
#' @import tibble
#' @import lubridate
#'
me2_read_residuals <- function(me2_res_file,
                               block_boundaries = list("start" = "^\\s*Scaled Residuals"),
                               species = NA,
                               dates = NA,
                               tz = "Etc/GMT-1") {

  # check if file me2_res_file
  if (!file.exists(me2_res_file)) {
    cli::cli_abort(c(
      "{.var me2_res_file} must be a valid file:",
      "x" = "File '{me2_res_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_res_file)
  text_filter <- text[text != ""]

  # residuals are located between "Scaled Residuals"
  index_start <- stringr::str_which(text_filter, block_boundaries$start)
  if (length(index_start) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The start header of the data block could not be found in '{me2_res_file}'.",
      "x" = "Did you provide the correct file?"
    ))
  }
  if (length(index_start) > 1) {
    index_end <- c(index_start[2:length(index_start)], length(text_filter) + 1) # need to add one line
  } else {
    index_end <- length(text_filter) + 1 # need to add one line
  }

  res_matrix <- tibble()
  for (run.number in seq(1, length(index_start), 1)) {
    matrix.text <- text_filter[seq(index_start[run.number] + 1, index_end[run.number] - 1, 1)]
    matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
      delim = " ",
      trim_ws = TRUE,
      col_names = FALSE
    )) %>%
      select(-dplyr::last(names(.))) %>%
      tibble::add_column(residual_type = "residual_scaled", .before = "X1") %>%
      tibble::add_column(base_run = run.number, .before = "X1")

    # we need to add the species names and datetime to the residuals.
    # set the names of the columns to the correct species
    if (!is.na(species)) {
      if (length(species) == (ncol(matrix.tmp) - 2)) {
        names(matrix.tmp)[3:ncol(matrix.tmp)] <- species
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the matrix with residuals ({ncol(matrix.tmp)-2})",
          "x" = "{.var species} must have the same length as the number of
          cols in the matrix containing residuals."
        ))
      }
    } else {
      # add default
      species_default <- paste0(
        "species_",
        sprintf(
          "%02d",
          seq(1, (ncol(matrix.tmp) - 2), 1)
        )
      )
      names(matrix.tmp)[3:ncol(matrix.tmp)] <- species_default
    }

    if (!is.na(dates)) {
      if (length(dates) == nrow(matrix.tmp)) {
        matrix.tmp <- matrix.tmp %>%
          tibble::add_column(
            identifier = dates,
            .after = "base_run"
          )
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var dates} has a different length ({length(dates)})
           compared to the matrix with residuals ({nrow(matrix.tmp)})",
          "x" = "{.var dates} must have the same length as the number of
          rows in the matrix containing the residuals"
        ))
      }
    } else {
      # calculate dates
      date_start <- lubridate::ymd("1970-01-01")
      date_end <- date_start + (nrow(matrix.tmp) - 1)
      dates_default <- lubridate::ymd(seq(
        from = date_start,
        to = date_end,
        by = 1
      ),
      tz = tz
      )
      # add dates column as identifier
      matrix.tmp <- matrix.tmp %>%
        tibble::add_column(
          date = dates_default,
          .after = "base_run"
        )
    }

    if (nrow(res_matrix) > 0) {
      res_matrix <- dplyr::bind_rows(res_matrix, matrix.tmp)
    } else {
      res_matrix <- matrix.tmp
    }
    rm(matrix.tmp)
  }

  return(res_matrix)
}
