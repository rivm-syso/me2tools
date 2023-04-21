#' Get the residuals from the residual file
#'
#' Typical use of Multilinear Engine version 2 (ME-2) will provide three
#' different files after a run: \dQuote{.dat}: machine readable format,
#' \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file with
#' auxiliary information (i.e. headers). This function reads the residual file
#' and supports two different formats. The first format is based on the EPA_PMF
#' format which consists of a block of residuals and a block of scaled residuals
#' for each run. The second format is based on the multi-time script, which only
#' provides the scaled residuals. In this file the blocks are separated by a
#' header row. The contents of this header row can be provided through the
#' \code{block_boundaries} parameter. This parameter is ignored in cae of 
#' reading the EPA-PMF formated results.
#' The column \code{residual_type} in the output denotes the type of residual,
#' being \dQuote{residual} and \dQuote{residual_scaled}.
#'
#' @param me2_res_file ME2 output file (.rsd), containing the results for the
#'   residuals.
#' @param tidy_output Should the output of both G and F be reshaped into tidy
#'   data? Default: FALSE
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
                               tidy_output = FALSE,
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
  
  # There are two types of residual file formats, one containing a header and
  # one containing two blocks of data (residuals and scaled residuals)
  # separated by empty lines. The format of the latter consists of residuals
  # first, followed by the scaled residuals from a single run.
  #
  # To determine the format, we check if the first line contains a header using
  # the block boundaries.
  if(stringr::str_detect(text_filter[[1]], block_boundaries$start)) {
    
    # residuals are located between "Scaled Residuals"
    index_start <- stringr::str_which(text_filter, block_boundaries$start)
    
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
        tibble::add_column(model_run = run.number, .before = "X1")
      
      # we need to add the species names and datetime to the residuals.
      # set the names of the columns to the correct species
      if (length(species) > 1) {
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
        if (is.na(species)) {
          # add default
          species_default <- paste0(
            "species_",
            sprintf(
              "%02d",
              seq(1, (ncol(matrix.tmp) - 2), 1)
            )
          )
          names(matrix.tmp)[3:ncol(matrix.tmp)] <- species_default
        } else {
          cli::cli_abort(c(
            "Different lengths:",
            "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({ncol(matrix.tmp) - 2})",
            "x" = "{.var species} must have the same length as the number of
            colums in the residual matrix or should be set to 'NA'"
          ))
        }
      }
      
      if (length(dates) > 1) {
        if (length(dates) == nrow(matrix.tmp)) {
          matrix.tmp <- matrix.tmp %>%
            tibble::add_column(
              date = dates,
              .after = "model_run"
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
        if (is.na(dates)) {
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
              .after = "model_run"
            )
        } else {
          cli::cli_abort(c(
            "Different lengths:",
            "i" = "{.var dates} has a different length ({length(dates)})
           compared to the matrix with residuals ({nrow(matrix.tmp)})",
            "x" = "{.var dates} must have the same length as the number of
          rows in the matrix containing the residuals or set as `NA`"
          ))
        }
      }
      
      if (nrow(res_matrix) > 0) {
        res_matrix <- dplyr::bind_rows(res_matrix, matrix.tmp)
      } else {
        res_matrix <- matrix.tmp
      }
      rm(matrix.tmp)
    }
  } else {
    # Technically the header block has not been detected so this might be
    # the second format with multiple blocks separated by empty lines.

    # the data block starts at 1, and then after each empty line
    start_blocks <- c(1, which(stringr::str_length(text) == 0) + 1)
    num_dates <- (start_blocks[[2]] - start_blocks[[1]]) - 1 # one empty line
    end_blocks <- start_blocks + num_dates - 1
    
    # due to additional empty lines at the end, blocks can be defined outside
    # the file limit.
    end_blocks <- end_blocks[which(end_blocks < length(text))]
    start_blocks <- start_blocks[seq(1, length(end_blocks), 1)]
    
    # There are two blocks at odd (residuals) and even (scaled residuals)
    res_indices <- seq(
      from = 1,
      to = length(start_blocks),
      by = 2
    )
    
    scaledres_indices <- seq(
      from = 2,
      to = length(start_blocks),
      by = 2
    )
    res_matrix <- tibble()
    
    # read the data
    for (run.number in seq(1, length(start_blocks)/2, 1)) {
      
      # get residual
      line_start <- start_blocks[[res_indices[[run.number]]]]
      line_end <- end_blocks[[res_indices[[run.number]]]]
      
      matrix.text <- text[seq(line_start, line_end, 1)]
      
      # transform to tibble
      matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
                                                             delim = " ",
                                                             trim_ws = TRUE,
                                                             col_names = FALSE
      ))  %>%
        select(-dplyr::first(names(.))) %>%
        tibble::add_column(residual_type = "residual", .before = "X2") %>%
        tibble::add_column(model_run = run.number, .before = "X2")
      
      
      # we need to add the species names and datetime to the residuals.
      # set the names of the columns to the correct species
      
      if (length(species) > 1) {
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
        if (is.na(species)) {
          # add default
          species_default <- paste0(
            "species_",
            sprintf(
              "%02d",
              seq(1, (ncol(matrix.tmp) - 2), 1)
            )
          )
          names(matrix.tmp)[3:ncol(matrix.tmp)] <- species_default
        } else {
          cli::cli_abort(c(
            "Different lengths:",
            "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({ncol(matrix.tmp) - 2})",
            "x" = "{.var species} must have the same length as the number of
            colums in the residual matrix or should be set to 'NA'"
          ))
        }
      }
      
      if (length(dates) > 1) {
        if (length(dates) == nrow(matrix.tmp)) {
          matrix.tmp <- matrix.tmp %>%
            tibble::add_column(
              date = dates,
              .after = "model_run"
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
        if (is.na(dates)) {
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
              .after = "model_run"
            )
        } else {
          cli::cli_abort(c(
            "Different lengths:",
            "i" = "{.var dates} has a different length ({length(dates)})
           compared to the matrix with residuals ({nrow(matrix.tmp)})",
            "x" = "{.var dates} must have the same length as the number of
          rows in the matrix containing the residuals or set as `NA`"
          ))
        }
      }
      
      # add data to output
      if (nrow(res_matrix) > 0) {
        res_matrix <- dplyr::bind_rows(res_matrix, matrix.tmp)
      } else {
        res_matrix <- matrix.tmp
      }
      rm(matrix.tmp)
      
      ########################################################
      # scaled residuals  
      ########################################################
      line_start <- start_blocks[[scaledres_indices[[run.number]]]]
      line_end <- end_blocks[[scaledres_indices[[run.number]]]]
      
      matrix.text <- text[seq(line_start, line_end, 1)]
      
      # transform to tibble
      matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
                                                             delim = " ",
                                                             trim_ws = TRUE,
                                                             col_names = FALSE
      ))  %>%
        select(-dplyr::first(names(.))) %>%
        tibble::add_column(residual_type = "residual_scaled", .before = "X2") %>%
        tibble::add_column(model_run = run.number, .before = "X2")
      
      
      # we need to add the species names and datetime to the residuals.
      # set the names of the columns to the correct species
      if (length(species) > 1) {
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
        if (is.na(species)) {
          # add default
          species_default <- paste0(
            "species_",
            sprintf(
              "%02d",
              seq(1, (ncol(matrix.tmp) - 2), 1)
            )
          )
          names(matrix.tmp)[3:ncol(matrix.tmp)] <- species_default
        } else {
          cli::cli_abort(c(
            "Different lengths:",
            "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({ncol(matrix.tmp) - 2})",
            "x" = "{.var species} must have the same length as the number of
            colums in the residual matrix or should be set to 'NA'"
          ))
        }
      }
      
      if (length(dates) > 1) {
        if (length(dates) == nrow(matrix.tmp)) {
          matrix.tmp <- matrix.tmp %>%
            tibble::add_column(
              date = dates,
              .after = "model_run"
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
        if (is.na(dates)) {
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
              .after = "model_run"
            )
        } else {
          cli::cli_abort(c(
            "Different lengths:",
            "i" = "{.var dates} has a different length ({length(dates)})
           compared to the matrix with residuals ({nrow(matrix.tmp)})",
            "x" = "{.var dates} must have the same length as the number of
          rows in the matrix containing the residuals or set as `NA`"
          ))
        }
      }
      # add data to output
      if (nrow(res_matrix) > 0) {
        res_matrix <- dplyr::bind_rows(res_matrix, matrix.tmp)
      } else {
        res_matrix <- matrix.tmp
      }
      rm(matrix.tmp)
    }
  }
  
  if (tidy_output) {
    species.order <- names(res_matrix[4:length(res_matrix)])
    
    res_matrix <- res_matrix %>% 
      pivot_longer(cols = -c("residual_type", "model_run", "date"),
                   names_to = "species",
                   values_to = "value") %>% 
      mutate(species = factor(species, levels = species.order))
  }

  return(res_matrix)
}
