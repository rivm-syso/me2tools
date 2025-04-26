#' Export G-matrix to MetCor
#'
#' This function exports the G-matrix to a file that can be used as input for
#' MetCor.
#'
#' @param G_matrix G-matrix that should be exported to a MetCor input file. This
#'   is typically the G_matrix provided by \code{\link{me2_read_G}} with the
#'   \code{tidy_output} parameter set to \code{FALSE}. This output should be
#'   expanded with the latitude (column: \code{lat}) and longitude (column:
#'   \code{lon}) coordinates of the receptor site(s). These coordinates should
#'   be exactly the same as the coordinates used in the **HYSPLIT** model for
#'   the receptor sites. Based on the unique lat and lon coordinates the
#'   number of receptor sites is automatically derived.
#' @param file File name and location to save the MetCor input file
#' @param time_res The time resolution between the measurements (in minutes),
#'   which is used to calculate the start (ITIME) and end times (FTIME) of each
#'   sample needed in MetCor. These times are used to associate the correct
#'   trajectories. If this variable is set to "auto" the function will try to
#'   derive the time resolution automatically from the data. Note: the "auto"
#'   feature only works when data sets of different time resolutions are not
#'   mixed.
#' @param time_format This parameter controls the definition of the start and
#'   end times in the MetCor output. \code{"0000"} ensures that the hour starts
#'   at '0000' and ends at '0059'. When using \code{"0001"} the hour starts at
#'   '0001' and ends at '0100'. Default is set to \code{"0000"}.
#' @param factor_names The names of the columns in the G-matrix that are
#'   identified as factors. These columns will be used to add the
#'   data/threshold values to output file. If \code{factor_names} is set to
#'   "auto" the function tries to find the factors based on the known columns
#'   in the G-matrix. Note, this feature will not work if non-default columns
#'   are present in the G-matrix.
#' @param na.rm remove rows with NA values for the factor contributions?
#'   Default is \code{TRUE}
#' @param neg.replace replace any negative value in the G contributions with 0.
#'   This is set to \code{TRUE} by default.
#' @param force_utc With this parameter the date time output for Metcor can be
#'   forced into UTC (the same timezone as Hysplit results). When set to
#'   \code{TRUE} there has to be no correction for the "Correlated Data Time
#'   Zone" within Metcor. If this parameter is \code{FALSE}, there has to be
#'   a manual correction when the time zone of the data is not UTC. Default
#'   setting of this parameter is \code{FALSE}, so no correction is applied when
#'   creating the output.
#'
#' @return tab-delimited plain text file containing date ranges, variable
#'   names, data/threshold values, receptor site coordinates and an optional
#'   multi site variable RECEPTORMAX if the number of receptor sites based on
#'   the uniqueness of the lat/lon coordinates in the input is more than 1.
#'
#' It should be noted that the full RECEPTORMAX, denoting the number of
#' receptor sites used in the analysis, is not required in MetCor. It is
#' perfectly fine to use a RECEPTORMAX smaller than the number of sites. In
#' this case, the user is advised to make manual changes to the output file.
#' The automatic default of RECEPTORMAX in MetCor is 2, when the RECEPTORMAX
#' is not given but the data is multi-site.
#'
#' @export
#'
#' @seealso \code{\link{metcor_import}}, \code{\link{metcor_project_raster}},
#' \code{\link{metcor_plot}}, \code{\link{metcor_plot_options}},
#' \code{\link{metcor_fix_hysplit}}
#'
#'
metcor_export <- function(G_matrix,
                          file,
                          time_res = "auto",
                          time_format = "0000",
                          factor_names = "auto",
                          na.rm = TRUE,
                          neg.replace = TRUE,
                          force_utc = FALSE) {

  # IDATE, ITIME, FDATE FTIME, LATR, LONR, FACTORS
  # LATR en LONR has to be exactly the same as the ones used in HySplit
  # IDATE and FDATE is the date in "yyyymmdd" format
  # ITIME and FTIME depend on the time resolution of the provided date times in
  # the G-matrix.

  # check for presence of lat/lon and give an error if they do not exists.
  if (!("lat" %in% colnames(G_matrix)) || !("lon" %in% colnames(G_matrix))) {
    cli::cli_abort(c(
      "{.var G_matrix} doesn't contain coordinates:",
      "x" = "Columns 'lat' and 'lon' should be available in the {.var G_matrix}."
    ))
  }

  # check for date
  if (!("date" %in% colnames(G_matrix))) {
    cli::cli_abort(c(
      "{.var G_matrix} doesn't have date time stamps:",
      "x" = "Column 'date' should be available in the {.var G_matrix}."
    ))
  }

  if(!("POSIXct" %in% class(G_matrix$date))) {
    cli::cli_abort(c(
      "The {.var date} in {.var G_matrix} is not POSIXct:",
      "x" = "Column 'date' should be a POSIXct date in the {.var G_matrix}."
    ))
  }

  # reset date to UTC
  if (force_utc) {
    check_tz <- lubridate::tz(G_matrix$date)
    if (check_tz != "UTC") {
      G_matrix$date <- lubridate::with_tz(G_matrix$date,
                                          tzone = "UTC")
    }
  }

  # check for factor
  if (("factor" %in% colnames(G_matrix))) {
    G_matrix <- G_matrix %>%
      pivot_wider(id_cols = c("model_type",
                              "unit",
                              "model_run",
                              "run_type",
                              "date",
                              "lat",
                              "lon"),
                  names_from = factor,
                  values_from = value)

    cli::cli_warn(c(
      "{.var G_matrix} was in the tidied format:",
      "x" = "An attempt has made to pivot the data to the wide format."
    ))
  }

  # get the resolution from the data if time_res = "auto"
  if (time_res == "auto") {
    diff_minutes <- as.numeric(difftime(utils::tail(G_matrix$date, -1),
      utils::head(G_matrix$date, -1),
      units = "mins"
    ))
    if (min(diff_minutes) == max(diff_minutes)) {
      # all in same resolution
      time_res <- min(diff_minutes)
    } else {
      # some date/time is missing, leading to larger time resolutions
      # use the median
      time_res <- stats::median(diff_minutes)
      cli::cli_warn(c(
        "Time resolution in {.var G_matrix} is not constant:",
        "i" = "Some samples might have been intentionally left out. The time
        resolution is set to {time_res} minutes, which corresponds to the
        median."
      ))
    }
  }

  # get the factor_names from the data if factor_names = "auto".
  if (factor_names == "auto") {
    # remove the default names from the G-matrix
    default_col_names <- c(
      "model_type",
      "unit",
      "model_run",
      "run_type",
      "date",
      "lat",
      "lon"
    )
    factor_names <- names(G_matrix %>%
      select(-all_of(default_col_names)))
  }

  # check if we need to replace any negative values with 0 from the columns
  # defined in variable factor_names
  if (neg.replace) {
    G_matrix <- G_matrix %>%
      dplyr::mutate(across(factor_names, ~if_else(.x < 0, 0, .x)))
  }

  if(time_format == "0000") {
    # ITIME starts at whole hour
    # FTIME ends one minute earlier to prevent overlap of time periods
    metcor.set <- G_matrix %>%
      mutate(
        IDATE = format(date, "%Y%m%d"),
        ITIME = format(date, "%H%M"),
        FDATE = format(date + lubridate::minutes(time_res - 1), "%Y%m%d"),
        FTIME = format(date + lubridate::minutes(time_res - 1), "%H%M")
      )
  } else if(time_format == "0001") {
    # ITIME starts one minute later to prevent overlap of time periods
    # FTIME ends on the complete time.
    metcor.set <- G_matrix %>%
      mutate(
        IDATE = format(date, "%Y%m%d"),
        ITIME = format(date + lubridate::minutes(1), "%H%M"),
        FDATE = format(date + lubridate::minutes(time_res), "%Y%m%d"),
        FTIME = format(date + lubridate::minutes(time_res), "%H%M")
      )
  } else {
    # unsupported time format
    cli::cli_abort(c(
      "The {.var time_format} as an unknown value:",
      "x" = "{.var time_format} should be either '0000' or '0001'."
    ))
  }

  # perform other actions
  metcor.set <- metcor.set %>%
    rename(
      LATR = lat,
      LONR = lon
    ) %>%
    select(IDATE, ITIME, FDATE, FTIME, LATR, LONR, all_of(factor_names))

  # check if we need to remove the rows with NA values
  if (na.rm) {
    metcor.set <- metcor.set %>%
    filter(complete.cases(.))
  }

  # note: when removing dates using EPA-PMF options there might be NA values
  # for the factor, we only use complete.cases to make sure. When dates are
  # removed from the input prior to the analysis, this situation will not
  # occur.

  utils::write.table(metcor.set,
    file = file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )

  # unique combinations of lat/lon point to number of receptors:
  receptors <- unique(G_matrix[, c("lat", "lon")])

  if (nrow(receptors) > 1) {
    # add max receptor line
    write(paste0("RECEPTORMAX\t", nrow(receptors)),
      file = file,
      append = TRUE
    )
  }
}
