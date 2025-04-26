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
    
    # new method, using the default edition of readr (2)
    matrix.text <- stringr::str_squish(matrix.text)
    matrix.tmp_2 <-readr::read_delim(I(matrix.text),
                                     delim = " ",
                                     trim_ws = TRUE,
                                     col_names = FALSE,
                                     show_col_types = FALSE
    )  %>%
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
    
    # new method, using the default edition of readr (2)
    matrix.text <- stringr::str_squish(matrix.text)
    matrix.tmp_2 <-readr::read_delim(I(matrix.text),
                                     delim = " ",
                                     trim_ws = TRUE,
                                     col_names = FALSE,
                                     show_col_types = FALSE
    )  %>%
      select(-dplyr::first(names(.))) # drop the first column with only numbers
    
    matrix.tmp_2 <- matrix.tmp_2 %>%
      purrr::set_names(colnames)
  }
  
  # check if new method is identical to old method
  check_output <- all.equal(matrix.tmp, 
                            matrix.tmp_2, 
                            trim.levels = FALSE, 
                            check.attributes = FALSE)
  
  if(!check_output) {
    cli::cli_abort("The old and new output from `readr` are not equal!")
  }

  # return old output
  return(matrix.tmp)
}
