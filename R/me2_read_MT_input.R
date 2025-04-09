#' Get the multi time data
#'
#' The data from the multi time data file can be read using this function. It
#' also performs some clean up, allowing the use of data in R.
#'
#' @section Data format:
#'   The format of the multi time data has to consists of the following columns,
#'   in the same order as given in the table.
#'
#'   \tabular{ll}{
#'   \strong{Name}    \tab \strong{Description}\cr
#'   Datestart:       \tab Start date of the measurement in the 'yyyy-mm-dd' or 'dd-mm-yyyy' format\cr
#'   Time             \tab Start time of the measurement in the 'hh:mm' format\cr
#'   Dateend:         \tab End date of the measurement in the 'yyyy-mm-dd' or 'dd-mm-yyyy' format\cr
#'   Time             \tab End time of the measurement in the 'hh:mm' format\cr
#'   Tzone            \tab Timezone associated with the start and end date times\cr
#'   Begin            \tab Unique number, expressing the beginning of the sample in time-units used in the data\cr
#'   Length           \tab Unique number, expressing the length of the sampling time in time-units used in the data\cr
#'   End              \tab Unique number, expressing the end of the sample in time-units used in the data\cr
#'   ...              \tab Columns of species with their appropriate names\cr
#'   }
#'
#'   For more information on the multi time data format see, for example,
#'   \doi{10.1016/j.scitotenv.2022.157981}
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
me2_read_MT_input <- function(file,
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
  input.data.org <- utils::read.table(
    file = file,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE
  )

  input.data <- input.data.org %>%
    rename(
      date.start.org = Datestart.,
      time.start.org = Time,
      date.end.org = Dateend.,
      time.end.org = Time.1,
      tzone = Tzone
    ) %>%
    select(-X0) # row number


  date_start <- lubridate::ymd(
    seq(
      from = lubridate::ymd("1970-01-01"),
      to = lubridate::ymd("1970-01-01") + (nrow(input.data) - 1),
      by = 1
    ),
    tz = tz
  )
  date_end <- date_start

  for (i in seq(1, nrow(input.data), by = 1)) {
    if (trimws(input.data$time.start.org[i]) == "") {
      input.data$time.start.org[i] <- "00:00"
    }
    if (!is.na(suppressWarnings(lubridate::ymd_hm(paste(trimws(input.data$date.start.org[i]), trimws(input.data$time.start.org[i])))))) {
      date_start[i] <- lubridate::ymd_hm(paste(trimws(input.data$date.start.org[i]), trimws(input.data$time.start.org[i])), tz = tz)
    } else if (!is.na(suppressWarnings(lubridate::dmy_hm(paste(trimws(input.data$date.start.org[i]), trimws(input.data$time.start.org[i])))))) {
      date_start[i] <- lubridate::dmy_hm(paste(trimws(input.data$date.start.org[i]), trimws(input.data$time.start.org[i])), tz = tz)
    } else if (!is.na(suppressWarnings(lubridate::mdy_hm(paste(trimws(input.data$date.start.org[i]), trimws(input.data$time.start.org[i])))))) {
      date_start[i] <- lubridate::mdy_hm(paste(trimws(input.data$date.start.org[i]), trimws(input.data$time.start.org[i])), tz = tz)
    } else {
      cli::cli_abort(c(
        "Datetime must be in 'yyyy-mm-dd hh:mm' or 'dd-mm-yyyy hh:mm' format:",
        "x" = "The datetime formats in '{file}' are not supported."
      ))
    }
    if (trimws(input.data$time.end.org[i]) == "") {
      input.data$time.end.org[i] <- "00:00"
    }
    if (!is.na(suppressWarnings(lubridate::ymd_hm(paste(trimws(input.data$date.end.org[i]), trimws(input.data$time.end.org[i])))))) {
      date_end[i] <- lubridate::ymd_hm(paste(trimws(input.data$date.end.org[i]), trimws(input.data$time.end.org[i])), tz = tz)
    } else if (!is.na(suppressWarnings(lubridate::dmy_hm(paste(trimws(input.data$date.end.org[i]), trimws(input.data$time.end.org[i])))))) {
      date_end[i] <- lubridate::dmy_hm(paste(trimws(input.data$date.end.org[i]), trimws(input.data$time.end.org[i])), tz = tz)
    } else if (!is.na(suppressWarnings(lubridate::mdy_hm(paste(trimws(input.data$date.end.org[i]), trimws(input.data$time.end.org[i])))))) {
      date_end[i] <- lubridate::mdy_hm(paste(trimws(input.data$date.end.org[i]), trimws(input.data$time.end.org[i])), tz = tz)
    } else {
      cli::cli_abort(c(
        "Datetime must be in 'yyyy-mm-dd hh:mm' or 'dd-mm-yyyy hh:mm' format:",
        "x" = "The datetime formats in '{file}' are not supported."
      ))
    }
  }

  input.data <- input.data %>%
    tibble::add_column(date = date_start, .before = "date.start.org") %>%
    tibble::add_column(date_end = date_end, .before = "date.start.org")

  # all in list
  input <- list("org" = input.data.org, "conc" = input.data.org, "unc" = input.data.org)
  input$conc <- input.data %>% select(
    -contains(unc_identifier)
  )
  input$unc <- input.data %>% select(
    date,
    date_end,
    date.start.org,
    time.start.org,
    date.end.org,
    time.end.org,
    tzone,
    Begin,
    Length,
    End,
    contains(unc_identifier)
  )
  # remove _std
  names(input$unc) <- names(input$conc)
  message("Please note that the uncertainties in this file might have been altered using the ME-2 script!")

  return(input)
}
