#' Get the minimum and maximum DISP values from the DISP results file
#'   (DISPres_x.txt)
#'
#' Files stored after a DISP run are named with a user-specific prefix, shown 
#' here as an asterisk (*). Three output files (*_DISP.dat,  *_DISP.txt and 
#' *_DISP.rsd) are stored after a DISP run. The function \code{me2_DISP_read_F} 
#' reads all the factor profiles in the text file and \code{me2_DISP_read_G} 
#' reads all the factor contributions. The residuals, stored in the 
#' *_DISP.rsd file, can be read using \code{ me2_read_residuals}. Besides these 
#' three files, four other files are produced as output, corresponding to 
#' dQmax = [4, 8, 16, 32] and are *_DISPres1.txt, *_DISPres2.txt, 
#' *_DISPres3.txt and *_DISPres4.txt. These files contain the minimum and 
#' maximum DISP results for a specific dQmax and any of these files can be read 
#' using this function. Each of these files have the same format, containing 
#' minimum and maximum DISP values in concentration units and as percentage of 
#' species units. These are all read and stored in the output.
#'
#' @param DISPres_file location and file name of the \dQuote{DISPres_?.txt} file,
#'   corresponding to either dQmax = [4, 8, 16, 32].
#' @param base_run The number of the base run associated with the DISP results
#' @param tidy_output Should the output be reshaped into tidy data? Default:
#'   FALSE
#' @param species A vector containing the names of the species for the rows in
#'   the F-matrix. If these species names are outputted in the ME-2 output as 
#'   the second column (a column of row numbers being the first), then these 
#'   values are used when \code{species = NA}. If this second column with names 
#'   is not available, all species are named as \dQuote{species_xx}, with xx 
#'   being an unique number starting at 1.
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
#' \preformatted{
#' mydata$factor <- factor(mydata$factor)
#' mydata$factor <- dplyr::recode_factor(mydata$factor,
#'                                       `factor_01` = "MyFirstName",
#'                                       `factor_02` = "MySecondName",
#'                                       ...
#' )
#' }
#'
#' Please note that the above will only work when the data is read with the
#' \code{tidy_output = TRUE} setting.
#'
#' @return \code{me2_DISP_read_res} returns an object of class ``me2tools''.
#'   The object includes five main components: \code{call}, the command used
#'   to read the data; \code{data}, the DISP data for each BS run;
#'   \code{F_format}, the aggregated DISP data in the same format as the
#'   F_matrix; \code{diag} the diagnostic line from the file, containing the error code (0
#'   is no error, 6 or 9 indicates that the run was aborted) and the largest
#'   drop of Q; and \code{swaps}, the swap counts for each dQmax level (4, 8, 
#'   15, 25). If retained, e.g., using 
#'   \code{output <- me2_BSDISP_read_res(file)}, this output can be used to
#'   recover the data, reproduce, or undertake further analysis.
#'
#'   An me2tools output can be manipulated using a number of generic operations,
#'   including \code{print}, \code{plot} and \code{summary}.
#'
#' @export
#' 
#' @noMd
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
me2_DISP_read_res <- function(DISPres_file,
                              base_run = 1,
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

  # add a progress bar
  cli::cli_progress_bar("Reading data", total = length(blocks))
  
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
        tibble::add_column(model_run = base_run, .before = "species")
    } else {
      matrix.tmp <- matrix.tmp %>%
        tibble::add_column(factor_profile = "percentage_of_species_sum", .before = "species") %>%
        tibble::add_column(model_run = base_run, .before = "species")
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
    
    # update progress bar
    cli::cli_progress_update()
  }
  
  # Set id variables
  id_variables <- c("factor_profile", "model_run", "species", "model_type", "run_type")
  
  # Make the table longer, so we can calculate the mid point
  matrix <- matrix %>%
    tidyr::pivot_longer(-dplyr::all_of(id_variables), 
                        names_to = "factor",
                        values_to = "value") %>%
    arrange(factor,
            factor_profile,
            species) %>% 
    tidyr::pivot_wider(id_cols = dplyr::all_of(c(id_variables[1:4], "factor")),
                       names_from = "run_type",
                       values_from = "value") %>% 
    dplyr::mutate(DISP_avg = (DISP_max+DISP_min)/2) %>% 
    tidyr::pivot_longer(-dplyr::all_of(c(id_variables[1:4], "factor")), 
                        names_to = "run_type", 
                        values_to = "value") %>%
    arrange(run_type,
            factor_profile,
            factor,
            species) 
  

  if(!tidy_output) {
    data <- matrix  %>% 
      tidyr::pivot_wider(id_cols = dplyr::all_of(id_variables), 
                         names_from = "factor", 
                         values_from = "value") %>%
      arrange(run_type,
              factor_profile,
              species) 
    
  } else {
    data <- matrix
  }
  
  # Read diagnostics and swaps
  diag.line <- readr::read_table(file = DISPres_file, 
                                 col_names = FALSE,
                                 skip = 1, # there is an empty line
                                 n_max = 1,
                                 show_col_types = FALSE) %>% 
    purrr::set_names(c("error_code",
                       "max_decrease_Q"))
  
  swaps <- readr::read_table(file = DISPres_file, 
                             col_names = FALSE, 
                             skip = 2,  # there is an empty line
                             n_max = 4,
                             show_col_types = FALSE)%>%
    purrr::set_names(colnames)

  output <- list("data" = data,
                 "F_format" = matrix,
                 "diag" = diag.line,
                 "swaps" = swaps,
                 call = match.call())
  class(output) <- "me2tools"
  
  # close progress bar
  cli::cli_progress_done()
  
  return(output)
}

############################################################################
############################################################################
###                                                                      ###
###                         DEPRECATED FUNCTIONS                         ###
###                                                                      ###
############################################################################
############################################################################


#' Get the minimum and maximum DISP values from the DISP results file
#'   (DISPres_x.txt)
#'
#' This function has been renamed to [me2_DISP_read_res()] and will be 
#' deprecated in the near-future.
#'
#' @param DISPres_file location and file name of the \dQuote{DISPres_?.txt} file,
#'   corresponding to either dQmax = \[4, 8, 16, 32\].
#' @param tidy_output Should the output be reshaped into tidy data? Default:
#'   FALSE
#' @param species A vector containing the names of the species for the rows in
#'   the F-matrix. If these species name are outputted in the ME-2 output as the
#'   second column (a column of row numbers being the first), then these values
#'   are used when \code{species = NA}. If this second column with names is not
#'   available all species are named as \dQuote{species_xx}, with xx being an
#'   unique number starting at 1.
#'   
#' @export
#' 
#' @seealso [me2_DISP_read_res()]
#'
me2_DISP_read_minmax <- function(DISPres_file,
                                 tidy_output = FALSE,
                                 species = NA) {
  
  cli::cli_warn(c(
    "DEPRECATED FUNCTION:",
    "i" = "This function is deprecated. Please use `me2_DISP_read_res` instead."
  ))
  
  output <- me2_DISP_read_res(DISPres_file = DISPres_file,
                              tidy_output = tidy_output,
                              species = species)
  return(output)
}