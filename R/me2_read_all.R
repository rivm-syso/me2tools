#' Function to read all available ME2 output.
#'
#' @param me2_txt_file ME2 output file (.txt), containing the results and
#'   auxiliary information.
#' @param me2_res_file ME2 output file (.rsd), containing the results for the
#'   residuals.
#' @param me2_log_file file name of the ME-2 log file.
#' @param species A vector containing the names of the species for the rows in
#'   the F-matrix and residual results. If these species name are outputted in
#'   the ME-2 output as the second column (a column of row numbers being the
#'   first), then these values are used when \code{species = NA}. If this second
#'   column with names is not available all species are named as
#'   \dQuote{species_xx}, with xx being an unique number starting at 1.
#' @param dates A vector containing the sample dates for the rows in the
#'   G-matrix and residual results. If these dates are outputted in the
#'   ME-2 output for the G-matrix as the second column (a column of row numbers
#'   being the first), then these values are used when \code{dates = NA}. The
#'   values inside the ME-2 output for the G-matrix are overwritten when
#'   \code{dates} is provided.
#' @param factor_mass This is a vector with the same length as the number of
#'   factor profiles. It contains the total mass in concentration units which is
#'   used to transform the G matrix from unity to concentration units.
#' @param rescale_unity In some cases the averages of the G factors are not
#'   equal to unity. By default a warning is given whenever this is the case.
#'   With this parameter set to \code{TRUE} each factor is multiplied by
#'   1/avg(factor), so that the G factors are scaled to unity again.
#' @param tidy_output Should the output of both G and F be reshaped into tidy
#'   data? Default: FALSE
#' @param dc_species character vector containing the species that should be
#'   subtracted from the missing mass to account for double counting. See also
#'   the section about double counting.
#' @param tz Parameter to control the timezone when parameter \dQuote{dates} is
#'   not used. Default: 'Etc/GMT-1'
#'
#' @return list with Q_values, F_matrix, G_matrix, Q_contrib, Residuals and LOG
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
#' @export
#'
#' @seealso \code{\link{me2_read_F}}, \code{\link{me2_read_G}},
#' \code{\link{me2_read_Q}}, \code{\link{me2_read_Qcontrib}},
#' \code{\link{me2_read_residuals}}, \code{\link{me2_read_log}},
#' \code{\link{me2_read_dat}}
#'

#'
me2_read_all <- function(
    me2_txt_file,
    me2_res_file,
    me2_log_file,
    species = NA,
    dates = NA,
    factor_mass = NA,
    rescale_unity = NA,
    tidy_output = FALSE,
    dc_species = NA,
    tz = "Etc/GMT-1") {

  # reading the files

  # Q
  Q_values <- me2_read_Q(me2_txt_file,
                         block_boundaries = list("start" = "^\\s*Sum-of-squares  Q, Qmain, Qaux"))

  Q_contrib <- me2_read_Qcontrib(me2_txt_file,
                                 block_boundaries = list(
                                   "start" = "^\\s*Contribution in Qmain by different variables",
                                   "end" = "^\\s*Task #,  seed:"
                                 ))
  # F
  F_matrix <- me2_read_F(me2_txt_file,
                         tidy_output = tidy_output,
                         block_boundaries = list("start" = "^\\s*Composition factors",
                                                 "end" = "^\\s*Contribution in Qmain by different variables"),
                         species = species,
                         dc_species = dc_species)


  # G
  G_matrix <- me2_read_G(me2_txt_file,
                         tidy_output = tidy_output,
                         block_boundaries = list("start" = "^\\s*Time factors",
                                                 "end" = "^\\s*Composition factors"),
                         dates = dates,
                         factor_mass = factor_mass,
                         rescale_unity = rescale_unity,
                         tz = tz)

  # Log file
  LOG <- me2_read_log(me2_log_file,
                      block_boundaries = list("start" = "^\\s*Starting to process with the same .INI file next task:"))


  # Residuals
  Residuals <- me2_read_residuals(me2_res_file,
                                  block_boundaries = list("start" = "^\\s*Scaled Residuals"),
                                  species = species,
                                  dates = dates,
                                  tz = tz)

  # return output
  output <- list("Q_values" = Q_values,
                 "F_matrix" = F_matrix,
                 "G_matrix" = G_matrix,
                 "Q_contrib" = Q_contrib,
                 "Residuals" = Residuals,
                 "LOG" = LOG)
  class(output) <- "me2tools"
  return(output)
}
