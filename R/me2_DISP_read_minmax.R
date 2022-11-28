#' Get the minimum and maximum DISP values from the DISP results file
#'   (DISPres_x.txt)
#'
#' The minimum and the maximum DISP values are stored in 4 different files,
#' corresponding to dQmax = 4, 8, 16, 32. Each of these files have the same
#' format, containing minimum and maximum DISP values in concentration units and
#' as percentage of species units. These are all read and stored in the output.
#'
#' @param DISPres_file location and file name of the \dQuote{DISPres_x} file,
#'   corresponding to dQmax = 4, 8, 16, 32.
#' @param tidy_output Should the output be reshaped into tidy data? Default:
#'   FALSE
#' @param species A vector containing the names of the species for the rows in
#'   the F-matrix. If these species name are outputted in the ME-2 output as the
#'   second column (a column of row numbers being the first), then these values
#'   are used when \code{species = NA}. If this second column with names is not
#'   available all species are named as \dQuote{species_xx}, with xx being an
#'   unique number starting at 1.
#'
#' @section Adding species to min/max:
#' By far the easiest way to add \dQuote{species} to the DISP min/max results is
#' to provide them as input parameters. The species names can probably be found
#' in the original data input used for ME-2 calculations.
#'
#' Other ways to replace the \dQuote{species} might be by using conditional
#' replacement of the default names or the use of solutions involving
#' \code{rep} to replicate elements of vectors.
#'
#' @section Renaming factor profiles:
#' As the labels of the factor profiles are unknown, this routine reads them as
#' \dQuote{factor_xx}, with xx being a unique number and outputted as a
#' character array. The easiest way to work with this data is by re-coding this
#' vector as \code{factor}. Then the order of the factor profiles and the
#' correct names can be easily set using the following code.
#'
#' ```R
#' mydata$factor <- factor(mydata$factor)
#' mydata$factor <- dplyr::recode_factor(mydata$factor,
#'                                       `factor_01` = "MyFirstName",
#'                                       `factor_02` = "MySecondName",
#'                                       ...
#' )
#' ```
#'
#' Please note that the above will only work when the data is read with the
#' \code{tidy_output = TRUE} setting.
#'
#' @return (tidied) tibble containing the minimum and maximum values for the
#'   Best-fit with original values for the DISP runs. The output of the
#'   min/max-values for the DISP result is identical to the output of the
#'   F-values from the base runs. The only difference is that the DISP results
#'   have different values in the \dQuote{run_type} column. In this case this
#'   column contains either \dQuote{DISP_min} or \dQuote{DISP_max}.
#'
#' @export
#'
#' @seealso \code{\link{me2_read_F}}, \code{\link{me2_BS_read_F}}, 
#' \code{\link{me2_DISP_read_F}}, \code{\link{me2_read_all}}, 
#' \code{\link{me2_read_dat}}
#'
#' @import cli
#' @import readr
#' @import stringr
#' @import purrr
#' @import tibble
#' @import dplyr
#'
me2_DISP_read_minmax <- function(DISPres_file,
                                 tidy_output = FALSE,
                                 species = NA) {


  # check if file exists
  if (!file.exists(DISPres_file)) {
    cli::cli_abort(c(
      "{.var DISPres_file} must be a valid file:",
      "x" = "File '{DISPres_file}' does not exists."
    ))
  }

  text <- readr::read_lines(DISPres_file)

  # In the file, there is a line with two numbers, followed by four lines of
  # data. In the first line, the first value is an error code: 0 means no error;
  # 6 or 9 indicates that the run was aborted. If this first value is non-zero,
  # the DISP analysis results are considered invalid. The second value is the
  # largest observed drop of Q during DISP.
  #
  # Below these diagnostics there are four blocks of data, where each column is
  # a factor and each row a species: (1) the profile matrix lower bound, in
  # concentration units; (2) the profile matrix upper bound, in concentration
  # units; (3) the profile matrix lower bound, in % species units; (4) the
  # profile matrix upper bound, in % species.

  # find the locations of the empty lines
  index_empty_lines <- which(stringr::str_length(text) == 0) + 1 # add one to skip empty row

  first_block_line <- 9
  # find the block_lines for the other blocks
  start_blocks <- index_empty_lines[which(index_empty_lines >= 9)]
  num_species <- (start_blocks[[2]] - start_blocks[[1]]) - 2
  end_blocks <- start_blocks + num_species - 1

  blocks <- seq(1, 4, 1)

  matrix <- tibble()
  for (block in blocks) {
    line_start <- start_blocks[[block]]
    line_end <- end_blocks[[block]]

    matrix.text <- text[seq(line_start, line_end, 1)]

    # transform to tibble
    matrix.tmp <- readr::with_edition(1, readr::read_delim(matrix.text,
      delim = " ",
      trim_ws = TRUE,
      col_names = FALSE
    ))

    colnames <- paste0("factor_", sprintf("%02d", seq(1, ncol(matrix.tmp), 1)))


    # create matrix
    matrix.tmp <- matrix.tmp %>%
      purrr::set_names(colnames)

    # check if we need to add species
    ## Check length of species against length of data!
    if (length(species) > 1) {
      if (length(species) == nrow(matrix.tmp)) {
        matrix.tmp <- matrix.tmp %>%
          tibble::add_column(
            identifier = species,
            .before = "factor_01"
          )
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the matrix ({nrow(matrix.tmp)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the matrix"
        ))
      }
    } else {
      if (!is.na(species)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the matrix ({nrow(matrix.tmp)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the matrix or should be set to 'NA'"
        ))
      }
    }

    ## add some default columns
    # default F has column identifier, but DISP has not
    if (!any(grepl("identifier", colnames(matrix.tmp)))) {
      matrix.tmp <- matrix.tmp %>%
        tibble::add_column(
          identifier = paste0(
            "species_",
            sprintf(
              "%02d",
              seq(1, nrow(matrix.tmp), 1)
            )
          ),
          .before = "factor_01"
        )
    }

    ## add PMFR identifier columns
    matrix.tmp <- matrix.tmp %>%
      mutate(identifier = factor(identifier, levels = identifier)) %>%
      dplyr::rename(species = identifier) %>%
      tibble::add_column(model_type = "ME-2", .before = "species")

    if (block < 3) {
      # block 1 & 2 are concentrations
      matrix.tmp <- matrix.tmp %>%
        tibble::add_column(factor_profile = "concentration_of_species", .before = "species") %>%
        tibble::add_column(model_run = 1, .before = "species")
    } else {
      matrix.tmp <- matrix.tmp %>%
        tibble::add_column(factor_profile = "percentage_of_species_sum", .before = "species") %>%
        tibble::add_column(model_run = 1, .before = "species")
    }

    if ((block == 1) || (block == 3)) {
      matrix.tmp <- matrix.tmp %>%
        tibble::add_column(run_type = "DISP_min", .before = "species")
    }
    if ((block == 2) || (block == 4)) {
      matrix.tmp <- matrix.tmp %>%
        tibble::add_column(run_type = "DISP_max", .before = "species")
    }

    if (nrow(matrix) > 0) {
      matrix <- dplyr::bind_rows(matrix, matrix.tmp)
    } else {
      matrix <- matrix.tmp
    }
  }

  if(tidy_output) {
    # Set id variables
    id_variables <- c("factor_profile", "model_run", "species", "model_type", "run_type")

    # Make the table longer
    matrix <- matrix %>%
      tidyr::pivot_longer(-dplyr::all_of(id_variables), names_to = "factor") %>%
      arrange(factor,
              factor_profile,
              species)

  }

  return(matrix)
}
