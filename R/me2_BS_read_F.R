#' Get the F values from the BS results file (.txt)
#'
#' Typical use of Bootstrap using Multilinear Engine version 2 (ME-2) will
#' provide three different files after a run: \dQuote{.dat}: machine readable
#' format, \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file
#' with auxiliary information (i.e. headers). The headers in the text (.txt)
#' file are used in me2tools to split the data into several blocks. In this
#' function the F values are read using the following header information,
#' denoting the \dQuote{start} and \dQuote{end} lines of the block containing
#' the required data.
#'
#' @param me2_bs_txt_file ME2 output file (.txt), containing the results and
#'   auxiliary information for the BS runs.
#' @param tidy_output Should the output be reshaped into tidy data? Default:
#'   FALSE
#' @param block_boundaries A list containing the \dQuote{start} and \dQuote{end}
#'   string used to identify the boundaries of the block containing the values
#'   for the F-matrix.
#' @param species A vector containing the names of the species for the rows in
#'   the F-matrix. If these species name are outputted in the ME-2 output as the
#'   second column (a column of row numbers being the first), then these values
#'   are used when \code{species = NA}. If this second column with names is not
#'   available all species are named as \dQuote{species_xx}, with xx being an
#'   unique number starting at 1.
#' @param dc_species character vector containing the species that should be
#' subtracted from the missing mass to account for double counting. See also the
#' section about double counting.
#'
#'
#' @section Defining the block with data:
#'
#' The \dQuote{.dat} file contains a lot of empty lines, which are stripped
#' immediately after reading the file. This is important, as this ensures that
#' each block of data is contained withing two identifying strings.
#'
#' The function uses these lines to determine the start and end line in which
#' the data lies. As such, the lines between start and end are read and parsed.
#'
#' The identifying strings are provided using the \code{block_boundaries}, which
#' is a list containing the \dQuote{start} and \dQuote{end} string. The string
#' can be preceded by \code{^\\s*} which tells the function that the string can
#' begin with an undefined number of spaces.
#'
#' @section Adding species to F:
#' By far the easiest way to add \dQuote{species} to F is to provide them as
#' input parameters. The species names can probably be found in the original
#' data input used for ME-2 calculations.
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
#' Please note that the above will only work when the data is read with the
#' \code{tidy_output = TRUE} setting.
#'
#' @section Double counting in missing mass variable:
#'
#' In some ME-2 applications a variable containing the missing mass is used.
#' This variable often contains the difference between the total mass (i.e.
#' PM10) and the sum of constituents measured in the total mass. Typically,
#' in factor analysis it is assumed that the sum of constituents is
#' approximately equal to the total mass. In other words, the constituents are
#' responsible for most, if not all, of the total mass.
#'
#' If this is not the case, there are non-measured constituents that contribute
#' to the total mass. As such, the missing mass (denoted as \dQuote{m.mass}
#' in me2tools) is calculated by subtracting the sum of the constituents from
#' the total mass and introduced in the factor analysis.
#'
#' However, in the case of multi-time (MT) factor analysis it might be
#' difficult to estimate the \dQuote{m.mass} due to the various time resolutions
#' in the data. Often only part of the available constituents can be summed and
#' subtracted from total mass (if they have the same sampling interval/time
#' resolution), leaving an overestimation of the \dQuote{m.mass}. After factor
#' analysis, this overestimation in the \dQuote{m.mass} still exists. Summing
#' the results of the factor analysis for each constituents, including
#' \dQuote{m.mass}, will therefore also be an overestimation due to double
#' counting of some masses. As such, comparison against the total mass will
#' fail.
#'
#' To prevent this overestimation the \dQuote{m.mass} variable should be
#' corrected for double counting. A viable way to do this in MT-factor analysis
#' is to subtract the sum of constituents that should have been used in the
#' initial calculation of \dQuote{m.mass} but could not be applied due to
#' different time resolutions. Since the F matrix from the factor analysis
#' is provided with unified time resolutions, the mass associated with these
#' constituents in the F matrix can now be subtracted to provide a better
#' estimate of \dQuote{m.mass}.
#'
#' The variable \dQuote{dc_species} is a character vector containing
#' constituents that should be subtracted from the \dQuote{m.mass} variable to
#' account for double counting.
#'
#' @return (tidied) tibble containing F-values for multiple BS runs. The output
#'   of the F-values for the BS results is identical to the output of the
#'   F-values from the base runs. The only difference is that the BS results
#'   have different values in the \dQuote{run_type} column. In this case this
#'   column contains information of the type of bootstrap.
#'
#' @export
#'
#' @seealso \code{\link{me2_read_F}}, \code{\link{me2_BS_read_F}},
#' \code{\link{me2_DISP_read_F}}, \code{\link{me2_DISP_read_minmax}},
#' \code{\link{me2_read_all}}, \code{\link{me2_read_dat}}
#'

#' @import cli
#' @import readr
#' @import stringr
#' @import tibble
#' @import dplyr
#'
me2_BS_read_F <- function(me2_bs_txt_file,
                          tidy_output = FALSE,
                          block_boundaries = list(
                            "start" = "^\\s*Factor matrix BB*",
                            "end" = "^\\s*Factor matrix CC"
                          ),
                          species = NA,
                          dc_species = NA) {

  # check if file exists
  if (!file.exists(me2_bs_txt_file)) {
    cli::cli_abort(c(
      "{.var me2_bs_txt_file} must be a valid file:",
      "x" = "File '{me2_bs_txt_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_bs_txt_file)
  text_filter <- text[text != ""]

  # factor profiles are located between "   Factor matrix BB" and
  # "   Factor matrix CC" (note the spaces)
  index_start <- stringr::str_which(text_filter, block_boundaries$start)
  if (length(index_start) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The start header of the data block could not be found in '{me2_bs_txt_file}'.",
      "x" = "Did you provide the correct file or the correct start string using
      {var. block_boundaries}?"
    ))
  }
  index_end <- stringr::str_which(text_filter, block_boundaries$end)
  if (length(index_end) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The end header of the data block could not be found in '{me2_bs_txt_file}'.",
      "x" = "Did you provide the correct file or the correct end string using
      {var. block_boundaries}?"
    ))
  }

  # the length of the index_start will provide how many runs were performed.
  # we need to get all runs
  f_matrix <- tibble::tibble()
  for (run.number in seq(1, length(index_start), 1)) {
    f_matrix.tmp <- extract_me2_matrix_between(
      text,
      index_start[run.number],
      index_end[run.number],
      headers = FALSE
    )

    # check if we need to add species
    ## Check length of species against length of data!
    if (length(species) > 1) {
      if (length(species) == nrow(f_matrix.tmp)) {
        if (!any(grepl("identifier", colnames(f_matrix.tmp)))) {
          f_matrix.tmp <- f_matrix.tmp %>%
            tibble::add_column(
              identifier = species,
              .before = "factor_01"
            )
        } else {
          f_matrix.tmp <- f_matrix.tmp %>%
            mutate(identifier == species)
        }
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(f_matrix)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix"
        ))
      }
    } else {
      if (!is.na(species)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(f_matrix)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix or should be set to 'NA'"
        ))
      }
    }

    f_matrix.tmp <- f_matrix.tmp %>%
      mutate(identifier = factor(identifier, levels = identifier)) %>%
      tidy_me2_factors(
        run_number = run.number,
        dc_species = dc_species,
        tidy_output = tidy_output
      )

    # add column with information
    F.type <- trimws(gsub("Factor matrix BB",
                          "",
                          text_filter[index_start[run.number]]))
    bootstrap.type <- "unknown"
    if (F.type == "Best-fit with original values") {
      bootstrap.type <- "BS_base"
    }
    if (F.type == "Bootstrapped values") {
      bootstrap.type <- "BS"
    }
    if (F.type == "Bootstrapped/pulled values") {
      bootstrap.type <- "BS_pulled"
    }

    f_matrix.tmp <- f_matrix.tmp %>%
      tibble::add_column(run_type = bootstrap.type, .before = "species")

    if (nrow(f_matrix) > 0) {
      f_matrix <- dplyr::bind_rows(f_matrix, f_matrix.tmp)
    } else {
      f_matrix <- f_matrix.tmp
    }
    rm(f_matrix.tmp)
  }
  return(f_matrix)
}
