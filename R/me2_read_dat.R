#' Read G and F matrices from ME-2 .dat file
#'
#' Typical use of Multilinear Engine version 2 (ME-2) will provide three
#' different files after a run: \dQuote{.dat}: machine readable format,
#' \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file with
#' auxiliary information (i.e. headers).
#'
#' This function reads the results stored in the \dQuote{.dat} file. This data
#' consists of both the numerical matrices of both G and F (also in this order)
#' for each of the defined number of runs.
#'
#' Additional information is needed to complement the G and F matrices. The
#' G matrix (contributions) are often associated with sample dates. These
#' sample dates can be given by the user using the \code{dates} variable. If
#' no dates are given, the function adds a daily date, starting from
#' \dQuote{1970-01-01} for every row in the G matrices.
#'
#' For the F matrices the species names are often needed, and can be provided
#' through the \code{species} variable. If no species are provided, default
#' species names (i.e. \dQuote{species_xx}) are provided, with xx being a
#' unique number starting from 1.
#'
#' @param me2_dat_file ME2 output file (.dat), containing the G and F matrices
#' @param factor_mass What mass from the F-matrix should be used to transform
#'   the G contributions from unity to concentration units. This parameter can
#'   be an integer (row number), the name of the species, or vector with the 
#'   same length of the number of factors. In the latter case, the first factor 
#'   is multiplied with the first item in this vector and so on.
#'   Default: \code{NA}, meaning no G matrix in concentration units are
#'   calculated.
#' @param dates This vector contains the sampling dates associated with the
#'   G-matrix. It should have the same length as the number of rows in G.
#'   Default: NA, meaning that the output is populated with daily dates starting
#'   from \dQuote{1970-01-01}.
#' @param species This vector contains the species associated with the F-matrix.
#'   It should have the same length as the number of rows in F. Default: NA,
#'   meaning that the output for F is populated as \dQuote{species_xx}, with xx
#'   being a unique number.
#' @param tidy_output Should the output of both G and F be reshaped into tidy
#'   data? Default: FALSE
#' @param rescale_unity In some cases the averages of the G factors are not
#'   equal to unity. By default a warning is given whenever this is the case.
#'   With this parameter set to \code{TRUE} each factor is multiplied by
#'   1/avg(factor), so that the G factors are scaled to unity again. Default:
#'   FALSE
#' @param threshold_unity The threshold which needs to be exceeded to provide
#'   a warning message. Default: 0.01
#' @param tz Parameter to control the time zone when parameter \dQuote{dates} is
#'   not used. Default: 'Etc/GMT-1'
#'
#' @return \code{me2_read_dat} returns an object of class ``me2tools''.
#'   The object includes three main components: \code{call}, the command used
#'   to generate the plot; \code{F_matrix}, the F matrices for each run;
#'   and \code{G_matrix}, the G matrices for each run. If retained, e.g., using
#'   \code{output <- me2_read_dat(file)}, this output can be used to recover
#'   the data, reproduce, or rework the original plot or undertake further
#'   analysis.
#'
#'   An me2tools output can be manipulated using a number of generic operations,
#'   including \code{print}, \code{plot} and \code{summary}.
#'
#' @section Adding dates and species to G and F:
#' By far the easiest way to add \dQuote{dates} to G and \dQuote{species} to F
#' is to provide them as input parameters. In both cases these data can
#' probably be found in the original data input used for ME-2 calculations.
#'
#' Other ways to replace the \dQuote{dates} and \dQuote{species} might be by
#' using conditional replacement of the default names or the use of solutions
#' involving \code{rep} to replicate elements of vectors.
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
#' @export
#'
#' @seealso \code{\link{me2_read_F}}, \code{\link{me2_read_G}},
#' \code{\link{me2_read_all}}
#'

#' @import tibble
#' @import cli
#' @import dplyr
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom utils read.table
#' 
me2_read_dat <- function (me2_dat_file,
                          factor_mass = NA,
                          dates = NA,
                          species = NA,
                          tidy_output = FALSE,
                          rescale_unity = FALSE,
                          threshold_unity = 0.01,
                          tz = "Etc/GMT-1"
) {


  #################################################################
  ##             Read file and split into list items             ##
  #################################################################

  # check if file exists
  if (!file.exists(me2_dat_file)) {
    cli::cli_abort(c(
      "{.var me2_dat_file} must be a valid file:",
      "x" = "File '{me2_dat_file}' does not exists."
    ))
  }

  filelines <- trimws(readLines(me2_dat_file))

  split.vec <- function(vec, sep = 0) {
    is.sep <- vec == sep
    split(vec[!is.sep], cumsum(is.sep)[!is.sep])
  }

  splitfilelines <- split.vec(filelines, sep = "")

  listoutput <- lapply(splitfilelines,
                       function(x) utils::read.table(textConnection(paste(x, collapse = "\n"))))
  
  # check and remove possible blocks with only one line of data
  drop.item <- c()
  for(list_index in 1:length(listoutput)) {
    if (nrow(listoutput[[list_index]]) == 1) {
      # drop this from the list...
      drop.item <-c(drop.item, list_index)
    }
  }
  listoutput[drop.item] <- NULL
  
  ##################################################################
  ##                 Creating necessary variables                 ##
  ##################################################################
  # The listoutput contains list items and can start with different numbers
  # To make sure we have the right data, we always use the index, not the named
  # number

  # add a progress bar
  cli::cli_progress_bar("Reading data", total = length(listoutput))
  
  num_of_runs <- length(listoutput)/2
  g_indices <- seq(from = 1, to = length(listoutput), by = 2)
  f_indices <- seq(from = 2, to = length(listoutput), by = 2)

  #################################################################
  ##                      Processing F data                      ##
  #################################################################
  F_matrix <- tibble()
  run_number <- 1
  for(f_index in f_indices) {
    # get the f values
    tmp_f_tibble <- listoutput[[f_index]]
    # add columnnames
    names(tmp_f_tibble) <- paste0(
      "factor_",
      sprintf(
        "%02d",
        seq(1, ncol(tmp_f_tibble), 1)
      )
    )

    ## if species is provided, do a cbind with the tmp_f_tibble and call this
    ## column "identifier", so it works well with the cleanup.

    # check if we need to add species
    ## Check length of species against length of data!
    if (length(species) > 1) {
      if (length(species) == nrow(tmp_f_tibble)) {
        if (!any(grepl("identifier", colnames(tmp_f_tibble)))) {
          tmp_f_tibble <- tmp_f_tibble %>%
            tibble::add_column(
              identifier = species,
              .before = "factor_01"
            )
        } else {
          tmp_f_tibble <- tmp_f_tibble %>%
            mutate(identifier == species)
        }
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_f_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix"
        ))
      }
    } else {
      if (!is.na(species)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(tmp_f_tibble)})",
          "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix or should be set to 'NA'"
        ))
      }
    }

    # cleanup
    tmp_f_tibble <- tidy_me2_factors(F_matrix = tmp_f_tibble,
                                     run_number = run_number,
                                     tidy_output = tidy_output) %>%
      tibble::add_column(run_type = "base_run", .before = "species")

    if (nrow(F_matrix) == 0) {
      F_matrix <- tmp_f_tibble
    } else {
      F_matrix <- dplyr::bind_rows(F_matrix, tmp_f_tibble)
    }
    run_number <- run_number+1
    # update progress bar
    cli::cli_progress_update()
  }

  #################################################################
  ##                      Processing G data                      ##
  #################################################################


  G_matrix <- tibble()
  run_number <- 1

  for(g_index in g_indices) {
    # get the g values
    tmp_g_tibble <- listoutput[[g_index]]
    # num_factors
    num_factors <- ncol(tmp_g_tibble)
    # add columnnames
    names(tmp_g_tibble) <- paste0(
      "factor_",
      sprintf(
        "%02d",
        seq(1, num_factors, 1)
      )
    )

    ## if date is provided, do a cbind with the tmp_g_tibble and call this
    ## column "identifier", so it works well with the cleanup.

    # add dates
    ## Check length of date against length of data!
    if (length(dates) > 1) {
      if (length(dates) == nrow(tmp_g_tibble)) {
        if (!any(grepl("identifier", colnames(tmp_g_tibble)))) {
          tmp_g_tibble <- tmp_g_tibble %>%
            tibble::add_column(
              identifier = dates,
              .before = "factor_01"
            )
        } else {
          tmp_g_tibble <- tmp_g_tibble %>%
            mutate(identifier == dates)
        }
        tmp_g_tibble <- tmp_g_tibble %>%
          mutate(identifier = format(identifier, "%Y-%m-%d %H:%M"))
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var dates} has a different length ({length(dates)})
           compared to the G-matrix ({nrow(tmp_g_tibble)})",
          "x" = "{.var dates} must have the same length as the number of
          rows in the G-matrix"
        ))
      }
    } else {
      if (!is.na(dates)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var dates} has a different length ({length(dates)})
           compared to the G-matrix ({nrow(tmp_g_tibble)})",
          "x" = "{.var dates} must have the same length as the number of
          rows in the G-matrix or should be set as 'NA'"
        ))
      }
    }

    if (length(factor_mass)==1) {
      if (!is.na(factor_mass)) {
        if (is.integer(factor_mass)) {
          factor_mass <- listoutput[[f_indices[run_number]]][factor_mass,]  
        } else if (is.character(factor_mass)) {
          # get the F_matrix for the current run_number
          current.F <- F_matrix %>% 
            filter(factor_profile == "concentration_of_species",
                   model_run == run_number,
                   species == factor_mass)
          factor_mass <- current.F$value
        } else {
          cli::cli_abort(c(
            "Variable should be integer or character:",
            "i" = "{.var factor_mass} with length {1} should be an integer or character.",
            "x" = "Type of length 1 should be integer or character."
          ))
        }
      } else {
        factor_mass <- NA
      }
    } else {
      if (length(factor_mass)==num_factors) {
        factor_mass <- factor_mass
      } else {
        cli::cli_abort(c(
          "Length should be 1 or number of factors:",
          "i" = "{.var factor_mass} has length {length(factor_mass)}.",
          "x" = "Length should be 1 or equal to {num_factors}"
        ))
      }
    }

    # cleanup
    tmp_g_tibble <- tidy_me2_contributions(G_matrix = tmp_g_tibble,
                                           run_number = run_number,
                                           factor_mass = factor_mass,
                                           rescale_unity = rescale_unity,
                                           threshold_unity = threshold_unity,
                                           tidy_output = tidy_output,
                                           tz = tz) %>%
      tibble::add_column(run_type = "base_runs", .after = "model_run")

    if (nrow(G_matrix) == 0) {
      G_matrix <- tmp_g_tibble
    } else {
      G_matrix <- dplyr::bind_rows(G_matrix, tmp_g_tibble)
    }
    run_number <- run_number+1
    # update progress bar
    cli::cli_progress_update()
  }

  output <- list("F_matrix" = F_matrix,
                 "G_matrix" = G_matrix,
                 "total_mass" = factor_mass,
                 call = match.call())
  class(output) <- "me2tools"
  
  # close progress bar
  cli::cli_progress_done()
  
  return(output)
}

