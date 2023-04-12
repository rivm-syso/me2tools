#' Get the input data
#'
#' The data from the input data file can be read using this function. It
#' also performs some clean up, allowing the use of data in R.
#'
#' @section Data format:
#'   The format of the multi time data has to consists of the following columns,
#'   in the same order as given in the table.
#'
#'   \tabular{ll}{
#'   \strong{Name}    \tab \strong{Description}\cr
#'   DateTime:        \tab Start date of the measurement in the 'yyyy-mm-dd hh:mm(:ss)' or 'dd-mm-yyyy hh:mm(:ss)' format\cr
#'   ...              \tab Columns of species with their appropriate names, alternating concentrations and uncertainties for each specie (i.e. specie1, unc_specie1, specie2, ...)\cr
#'   }
#'
#' @param file File name of the multi time data file.
#' @param unc_identifier Vector containing the prefix/suffix to denote the
#'   uncertainty for a species. Examples are \code{xx_std}, \code{xx_unc},
#'   \code{U_xx} and so on. The function splits the data into a data and
#'   uncertainty set by splitting the columns that contain one or more of the
#'   \code{unc_identifier}. Default is
#'   \code{c(".std", "_std", "_unc", ".unc", "U_", "U.", "unc_", "unc.")}.
#' @param sep The column separator used to separate values on one line. The
#'   default is \code{"\t"} which is a TAB character. Other options are
#'   \code{";"} or \code{","}
#' @param tz Time zone for date time, defaults to \code{Etc/GMT-1}
#'
#' @return list containing full data, both concentrations and uncertainties.
#'
#' @export
#'
#' @import cli
#' @import lubridate
#' @import tibble
#' @importFrom utils read.table
#'
me2_read_input <- function(file,
                           unc_identifier = c(
                             ".std",
                             "_std",
                             "_unc",
                             ".unc",
                             "U_",
                             "U.",
                             "unc_",
                             "unc."
                           ),
                           sep = "\t",
                           tz = "Etc/GMT-1") {
  # check if file exists
  if (!file.exists(file)) {
    cli::cli_abort(c(
      "{.var file} must be a valid file:",
      "x" = "File '{file}' does not exists."
    ))
  }

  # we need to add the species names and datetime to the residuals.
  # This information can be extracted from the input file.
  input.data <- utils::read.table(
    file = file,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE
  )
  # always use default date_time name for first column.
  names(input.data)[1] <- "date_time"

  if (!is.na(suppressWarnings(lubridate::ymd_hm(input.data[1, 1])))) {
    input.data$date_time <- lubridate::ymd_hm(input.data$date_time, tz = tz)
  } else if (!is.na(suppressWarnings(lubridate::dmy_hm(input.data[1, 1])))) {
    input.data$date_time <- lubridate::dmy_hm(input.data$date_time, tz = tz)
  } else if (!is.na(suppressWarnings(lubridate::mdy_hm(input.data[1, 1])))) {
    input.data$date_time <- lubridate::mdy_hm(input.data$date_time, tz = tz)
  } else if (!is.na(suppressWarnings(lubridate::ymd_hms(input.data[1, 1])))) {
    input.data$date_time <- lubridate::ymd_hms(input.data$date_time, tz = tz)
  } else if (!is.na(suppressWarnings(lubridate::dmy_hms(input.data[1, 1])))) {
    input.data$date_time <- lubridate::dmy_hms(input.data$date_time, tz = tz)
  } else if (!is.na(suppressWarnings(lubridate::mdy_hms(input.data[1, 1])))) {
    input.data$date_time <- lubridate::mdy_hms(input.data$date_time, tz = tz)
  } else {
    cli::cli_abort(c(
      "Datetime must be in 'yyyy-mm-dd hh:mm' or 'dd-mm-yyyy hh:mm' format:",
      "x" = "The datetime formats in '{file}' are not supported."
    ))
  }

  # all in list
  input <- list("org" = input.data, "conc" = input.data, "unc" = input.data)
  class(input) <- "me2tools"
  # replace with actual values.
  input$conc <- input.data %>% select(
    -contains(unc_identifier)
  )
  input$unc <- input.data %>% select(
    date_time,
    contains(unc_identifier)
  )
  # remove unc_identifier
  names(input$unc) <- names(input$conc)
  message("Please note that the uncertainties in this file might have been altered using the ME-2 script!")

  return(input)
}
