#' Get the F values from the BSDISP results file (BSDISPres?.txt)
#'
#' Files stored after a BSDISP run are named with a user-specific prefix, shown
#' here as an asterisk (*). Three output files (*_BSDISP.dat,  *_BSDISP.txt and
#' *_BSDISP.rsd) are stored after the DISP phase of the BSDISP run (note that
#' the BS phase of the BS-DISP has outputs with a \dQuote{_BS_} suffix, i.e.,
#' *_BS_.dat, *_BS_.txt, and *_BS_.rsd) and currently there are no functions
#' available for reading information from the \dQuote{txt} file. The residuals,
#' stored in the *_BSDISP.rsd file, can be read using \code{me2_read_residuals}.
#' Besides these three files, four other files are produced as output,
#' corresponding to dQmax = [0.5, 1, 2, 4] and are: *_BSDISPres1.txt,
#' *_BSDISPres2.txt, *_BSDISPres3.txt and *_DISPres4.txt, which can be read with
#' this function. These files contain the results for each BS resample after
#' displacing both up and down. In this function, the results are also
#' aggregated to \dQuote{BSDISP_P05}, \dQuote{BSDISP_P95} and
#' \dQuote{BSDISP_avg}, based on the BSDISP results in concentration units.
#'
#' @param BSDISPres_file location and file name of the \dQuote{DISPres_?.txt} file,
#'   corresponding to either dQmax = [0.5, 1, 2, 4].
#' @param base_run The number of the base run associated with the BSDISP results
#' @param tidy_output Should the output be reshaped into tidy data? Default:
#'   FALSE
#' @param species A vector containing the names of the species for the rows in
#'   the F-matrix. If these species name are outputted in the ME-2 output as the
#'   second column (a column of row numbers being the first), then these values
#'   are used when \code{species = NA}. If this second column with names is not
#'   available all species are named as \dQuote{species_xx}, with xx being an
#'   unique number starting at 1.
#' @param dc_species character vector containing the species that should be
#'   subtracted from the missing mass to account for double counting. See also
#'   the section about double counting. This will only be applied to the
#'   concentration values
#'
#' @section Adding species names to F:
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
#' \preformatted{
#' mydata$factor <- factor(mydata$factor)
#' mydata$factor <- dplyr::recode_factor(mydata$factor,
#'                                       `factor_01` = "MyFirstName",
#'                                       `factor_02` = "MySecondName",
#'                                       ...
#' )
#' }
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
#' in the input data. Often only part of the available constituents can be summed 
#' and subtracted from total mass (if they have the same sampling interval/time
#' resolution), leading to an overestimation of the \dQuote{m.mass}. After factor
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
#' different time resolutions (e.g., sum of metals or any other species collected
#' on a lower time resolution)' from the m.mass value for that factor. Since the 
#' F matrix from the factor analysis is provided with unified time resolutions, 
#' the mass associated with these constituents in the F matrix can now be 
#' subtracted to provide a better estimate of \dQuote{m.mass}.
#'
#' The variable \dQuote{dc_species} is a character vector containing
#' constituents that should be subtracted from the \dQuote{m.mass} variable to
#' account for double counting.
#'
#' @return \code{me2_BSDISP_read_res} returns an object of class ``me2tools''.
#'   The object includes three main components: \code{call}, the command used
#'   to read the data; \code{data}, the BSDISP data for each BS run;
#'   and \code{F_format}, the aggregated BSDISP data in the same format as the
#'   F_matrix. If retained, e.g., using 
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
#' \code{\link{me2_DISP_read_F}}, \code{\link{me2_DISP_read_res}},
#' \code{\link{me2_read_all}}, \code{\link{me2_read_dat}}
#'

#' @import cli
#' @import readr
#' @import stringr
#' @import tibble
#' @import dplyr
#'
me2_BSDISP_read_res <- function(BSDISPres_file,
                                base_run = 1,
                                tidy_output = FALSE,
                                species = NA,
                                dc_species = NA) {
  #################################################################
  ##             Read file and split into list items             ##
  #################################################################

  # check if file exists
  if (!file.exists(BSDISPres_file)) {
    cli::cli_abort(c(
      "{.var BSDISPres_file} must be a valid file:",
      "x" = "File '{BSDISPres_file}' does not exists."
    ))
  }

  filelines <- trimws(readr::read_lines(
    file = BSDISPres_file,
    skip = 6
  ))

  split.vec <- function(vec, sep = 0) {
    is.sep <- vec == sep
    split(vec[!is.sep], cumsum(is.sep)[!is.sep])
  }

  splitfilelines <- split.vec(filelines, sep = "")

  listoutput <- lapply(
    splitfilelines,
    function(x) utils::read.table(textConnection(paste(x, collapse = "\n")))
  )

  ##################################################################
  ##                 Creating necessary variables                 ##
  ##################################################################
  # The listoutput contains list items and can start with different numbers
  # To make sure we have the right data, we always use the index, not the named
  # number

  num_of_runs <- length(listoutput) / 5

  # define the four different blocks of data
  disp_d_c_indices <- seq(
    from = 1,
    to = length(listoutput) - num_of_runs,
    by = 4
  )
  disp_u_c_indices <- seq(
    from = 2,
    to = length(listoutput) - num_of_runs,
    by = 4
  )
  disp_d_p_indices <- seq(
    from = 3,
    to = length(listoutput) - num_of_runs,
    by = 4
  )
  disp_u_p_indices <- seq(
    from = 4,
    to = length(listoutput) - num_of_runs,
    by = 4
  )


  #################################################################
  ##                 Processing DISP-D_conc data                 ##
  #################################################################

  DISP_D_conc <- tibble()
  disp_d_c_index <- disp_d_c_indices[1]

  # add a progress bar
  cli::cli_progress_bar("Reading data", total = length(listoutput) - num_of_runs)

  for (disp_d_c_index in disp_d_c_indices) {
    # get the f values
    tmp_tibble <- listoutput[[disp_d_c_index]] %>%
      dplyr::filter(row_number() <= n() - 1) # drop the last row

    # add columnnames
    names(tmp_tibble) <- c(
      "run_number_final",
      "dummy",
      paste0(
        "factor_",
        sprintf(
          "%02d",
          seq(1, ncol(tmp_tibble) - 2, 1)
        )
      )
    )

    ## if species is provided, do a cbind with the tmp_f_tibble and call this
    ## column "identifier", so it works well with the cleanup.

    # check if we need to add species
    ## Check length of species against length of data!
    if (length(species) > 1) {
      if (length(species) == nrow(tmp_tibble)) {
        if (!any(grepl("identifier", colnames(tmp_tibble)))) {
          tmp_tibble <- tmp_tibble %>%
            tibble::add_column(
              identifier = species,
              .before = "factor_01"
            )
        } else {
          tmp_tibble <- tmp_tibble %>%
            mutate(identifier == species)
        }
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix"
        ))
      }
    } else {
      if (!is.na(species)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix or should be set to 'NA'"
        ))
      }
    }

    # cleanup
    tmp_tibble <- tidy_me2_factors(
      F_matrix = tmp_tibble %>% select(-dummy),
      run_number = -999,
      dc_species = dc_species
    ) %>%
      dplyr::filter(factor_profile == "concentration_of_species") %>%
      tibble::add_column(run_type = "BSDISP-Down", .before = "species") %>%
      dplyr::mutate(
        model_run = run_number_final,
        factor_profile = "concentration_of_species"
      ) %>%
      select(-run_number_final)

    if (nrow(DISP_D_conc) == 0) {
      DISP_D_conc <- tmp_tibble
    } else {
      DISP_D_conc <- dplyr::bind_rows(DISP_D_conc, tmp_tibble)
    }

    # update progress bar
    cli::cli_progress_update()
  }

  #################################################################
  ##                 Processing DISP-U_conc data                 ##
  #################################################################

  DISP_U_conc <- tibble()
  disp_u_c_index <- disp_u_c_indices[1]

  for (disp_u_c_index in disp_u_c_indices) {
    # get the f values
    tmp_tibble <- listoutput[[disp_u_c_index]] %>%
      dplyr::filter(row_number() <= n() - 1) # drop the last row

    # add columnnames
    names(tmp_tibble) <- c(
      "run_number_final",
      "dummy",
      paste0(
        "factor_",
        sprintf(
          "%02d",
          seq(1, ncol(tmp_tibble) - 2, 1)
        )
      )
    )

    ## if species is provided, do a cbind with the tmp_f_tibble and call this
    ## column "identifier", so it works well with the cleanup.

    # check if we need to add species
    ## Check length of species against length of data!
    if (length(species) > 1) {
      if (length(species) == nrow(tmp_tibble)) {
        if (!any(grepl("identifier", colnames(tmp_tibble)))) {
          tmp_tibble <- tmp_tibble %>%
            tibble::add_column(
              identifier = species,
              .before = "factor_01"
            )
        } else {
          tmp_tibble <- tmp_tibble %>%
            mutate(identifier == species)
        }
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix"
        ))
      }
    } else {
      if (!is.na(species)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix or should be set to 'NA'"
        ))
      }
    }

    # cleanup
    tmp_tibble <- tidy_me2_factors(
      F_matrix = tmp_tibble %>% select(-dummy),
      run_number = -999,
      dc_species = dc_species
    ) %>%
      dplyr::filter(factor_profile == "concentration_of_species") %>%
      tibble::add_column(run_type = "BSDISP-Up", .before = "species") %>%
      dplyr::mutate(
        model_run = run_number_final,
        factor_profile = "concentration_of_species"
      ) %>%
      select(-run_number_final)

    if (nrow(DISP_U_conc) == 0) {
      DISP_U_conc <- tmp_tibble
    } else {
      DISP_U_conc <- dplyr::bind_rows(DISP_U_conc, tmp_tibble)
    }

    # update progress bar
    cli::cli_progress_update()
  }

  #################################################################
  ##                 Processing DISP-D_perc data                 ##
  #################################################################

  DISP_D_perc <- tibble()
  disp_d_p_index <- disp_d_p_indices[1]

  for (disp_d_p_index in disp_d_p_indices) {
    # get the f values
    tmp_tibble <- listoutput[[disp_d_p_index]] %>%
      dplyr::filter(row_number() <= n() - 1) # drop the last row

    # add columnnames
    names(tmp_tibble) <- c(
      "run_number_final",
      "dummy",
      paste0(
        "factor_",
        sprintf(
          "%02d",
          seq(1, ncol(tmp_tibble) - 2, 1)
        )
      )
    )

    ## if species is provided, do a cbind with the tmp_f_tibble and call this
    ## column "identifier", so it works well with the cleanup.

    # check if we need to add species
    ## Check length of species against length of data!
    if (length(species) > 1) {
      if (length(species) == nrow(tmp_tibble)) {
        if (!any(grepl("identifier", colnames(tmp_tibble)))) {
          tmp_tibble <- tmp_tibble %>%
            tibble::add_column(
              identifier = species,
              .before = "factor_01"
            )
        } else {
          tmp_tibble <- tmp_tibble %>%
            mutate(identifier == species)
        }
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix"
        ))
      }
    } else {
      if (!is.na(species)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix or should be set to 'NA'"
        ))
      }
    }

    # cleanup
    tmp_tibble <- tidy_me2_factors(
      F_matrix = tmp_tibble %>% select(-dummy),
      run_number = -999
    ) %>%
      dplyr::filter(factor_profile == "concentration_of_species") %>%
      tibble::add_column(run_type = "BSDISP-Down", .before = "species") %>%
      dplyr::mutate(
        model_run = run_number_final,
        factor_profile = "percentage_of_species"
      ) %>%
      select(-run_number_final)

    if (nrow(DISP_D_perc) == 0) {
      DISP_D_perc <- tmp_tibble
    } else {
      DISP_D_perc <- dplyr::bind_rows(DISP_D_perc, tmp_tibble)
    }

    # update progress bar
    cli::cli_progress_update()
  }

  #################################################################
  ##                 Processing DISP-U_perc data                 ##
  #################################################################

  DISP_U_perc <- tibble()
  disp_u_p_index <- disp_u_p_indices[1]

  for (disp_u_p_index in disp_u_p_indices) {
    # get the f values
    tmp_tibble <- listoutput[[disp_u_p_index]] %>%
      dplyr::filter(row_number() <= n() - 1) # drop the last row

    # add columnnames
    names(tmp_tibble) <- c(
      "run_number_final",
      "dummy",
      paste0(
        "factor_",
        sprintf(
          "%02d",
          seq(1, ncol(tmp_tibble) - 2, 1)
        )
      )
    )

    ## if species is provided, do a cbind with the tmp_f_tibble and call this
    ## column "identifier", so it works well with the cleanup.

    # check if we need to add species
    ## Check length of species against length of data!
    if (length(species) > 1) {
      if (length(species) == nrow(tmp_tibble)) {
        if (!any(grepl("identifier", colnames(tmp_tibble)))) {
          tmp_tibble <- tmp_tibble %>%
            tibble::add_column(
              identifier = species,
              .before = "factor_01"
            )
        } else {
          tmp_tibble <- tmp_tibble %>%
            mutate(identifier == species)
        }
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix"
        ))
      }
    } else {
      if (!is.na(species)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix or should be set to 'NA'"
        ))
      }
    }

    # cleanup
    tmp_tibble <- tidy_me2_factors(
      F_matrix = tmp_tibble %>% select(-dummy),
      run_number = -999
    ) %>%
      dplyr::filter(factor_profile == "concentration_of_species") %>%
      tibble::add_column(run_type = "BSDISP-Up", .before = "species") %>%
      dplyr::mutate(
        model_run = run_number_final,
        factor_profile = "percentage_of_species"
      ) %>%
      select(-run_number_final)

    if (nrow(DISP_U_perc) == 0) {
      DISP_U_perc <- tmp_tibble
    } else {
      DISP_U_perc <- dplyr::bind_rows(DISP_U_perc, tmp_tibble)
    }

    # update progress bar
    cli::cli_progress_update()
  }


  # combine all the data
  bs_disp_data <- dplyr::bind_rows(
    DISP_D_conc,
    DISP_U_conc,
    DISP_D_perc,
    DISP_U_perc
  )

  # Set id variables
  id_variables <- c("factor_profile", "model_run", "species", "model_type", "run_type")

  # calculate stats
  bs_disp_data_stats <- bs_disp_data %>%
    filter(factor_profile == "concentration_of_species") %>%
    pivot_longer(
      cols = -dplyr::all_of(id_variables),
      names_to = "factor",
      values_to = "value"
    ) %>%
    group_by(factor, species) %>%
    summarise(
      BSDISP_avg = mean(value),
      BSDISP_P05 = epa_percentile(value, prob = 0.05),
      BSDISP_P95 = epa_percentile(value, prob = 0.95)
    ) %>%
    ungroup(factor, species) %>%
    pivot_longer(
      cols = c("BSDISP_avg", "BSDISP_P05", "BSDISP_P95"),
      names_to = "run_type",
      values_to = "value"
    ) %>%
    select(run_type, species, factor, value) %>%
    tibble::add_column(model_type = "ME-2", .before = "run_type") %>%
    tibble::add_column(factor_profile = "concentration_of_species", .before = "run_type") %>%
    tibble::add_column(model_run = base_run, .before = "run_type")

  # tidy output
  if (tidy_output) {
    id_variables <- c("factor_profile", "model_run", "species", "model_type", "run_type")
    bs_disp_data <- bs_disp_data %>%
      tidyr::pivot_longer(
        cols = -dplyr::all_of(id_variables),
        names_to = "factor",
        values_to = "value"
      )
  }

  output <- list(
    "data" = bs_disp_data,
    "F_format" = bs_disp_data_stats,
    call = match.call()
  )
  class(output) <- "me2tools"

  # close progress bar
  cli::cli_progress_done()

  return(output)
}
