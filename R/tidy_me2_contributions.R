#' Tidy ME-2 contributions
#'
#' This function tidies the G-matrix to match the PMFR output. This is an
#' internal function.
#'
#' @param G_matrix The preliminary output from functions reading G values
#' @param factor_mass This is a vector with the same length as the number of
#'   factor profiles. It contains the total mass in concentration units which is
#'   used to transform the G matrix from unity to concentration units.
#' @param run_number The number associated with the ME-2 run. Default: 1
#' @param tidy_output Should the factor contributions be reshaped into tidy
#'   data? Default: FALSE
#' @param rescale_unity In some cases the average of the G factors are not equal
#'   to unity. By default a warning is given whenever this is the case.
#'   With this parameter set to \code{TRUE} each factor is multiplied by
#'   1/avg(factor), so that the G factors are scaled to unity again.
#' @param threshold_unity The threshold which needs to be exceeded to provide
#'   a warning message. Default: 0.01
#' @param tz Parameter to control the timezone when parameter \dQuote{dates} is
#'   not used. Default: 'Etc/GMT-1'
#'
#' @noRd
#'
#' @import dplyr
#' @import tibble
#' @import lubridate
#' @import cli
#'
tidy_me2_contributions <- function(G_matrix,
                                   factor_mass = NA,
                                   run_number = 1,
                                   tidy_output = FALSE,
                                   rescale_unity = FALSE,
                                   threshold_unity = 0.01,
                                   tz = "Etc/GMT-1") {

  ## check if there is a identifier column, if not, add date
  # default F has column identifier, but bootstrap has not
  if (!any(grepl("identifier", colnames(G_matrix)))) {
    # calculate dates
    date_start <- lubridate::ymd("1970-01-01")
    date_end <- date_start + (nrow(G_matrix) - 1)
    dates <- lubridate::ymd(seq(from=date_start, to=date_end, by = 1 ), tz = tz)
    # add dates column as identifier
    G_matrix <- G_matrix %>%
      tibble::add_column(date = dates,
        .before = "factor_01"
      )
  } else {
    # try to process default identifier column as date
    ## add PMFR identifier columns
    G_matrix <- G_matrix %>%
      dplyr::rename(date = identifier)

    if (!is.na(lubridate::ymd_hm(G_matrix$date[1]))) {
      G_matrix <- G_matrix %>% mutate(date = lubridate::ymd_hm(date, tz = tz))
    } else if (!is.na(lubridate::ymd_hms(G_matrix$date[1]))) {
      G_matrix <- G_matrix %>% mutate(date = lubridate::ymd_hms(date, tz = tz))
    } else if (!is.na(lubridate::dmy_hm(G_matrix$date[1]))) {
      G_matrix <- G_matrix %>% mutate(date = lubridate::dmy_hm(date, tz = tz))
    } else if (!is.na(lubridate::dmy_hms(G_matrix$date[1]))) {
      G_matrix <- G_matrix %>% mutate(date = lubridate::dmy_hms(date, tz = tz))
    } else if (!is.na(lubridate::ymd(G_matrix$date[1]))) {
      G_matrix <- G_matrix %>% mutate(date = lubridate::ymd(date, tz = tz))
    } else if (!is.na(lubridate::dmy(G_matrix$date[1]))) {
      G_matrix <- G_matrix %>% mutate(date = lubridate::dmy(date, tz = tz))
    } else {
      cli::cli_abort(c(
        "Unknown date format:",
        "i" = "The date format in {.var G_matrix} is not supported.",
        "x" = "Date format should be 'ymd'/'dmy' or 'ymd_hm(s)'/'dmy_hm(s)'."
      ))
    }
  }

  # add other default columns
  G_matrix <- G_matrix %>%
    tibble::add_column(model_type = "ME-2", .before = "date") %>%
    tibble::add_column(unit = "normalised", .before = "date") %>%
    tibble::add_column(model_run = run_number, .before = "date")


  # rescale G to unity again before applying total mass and use that as normalised
  if (identical(rescale_unity, TRUE)) {
    # get the current normalised values
    scaled.data <- G_matrix %>%
      filter(unit == "normalised") %>%
      select(contains("factor"))
    # rename the original normalised values
    G_matrix <- G_matrix %>%
      mutate(unit = ifelse(unit == "normalised", "normalised_original", unit))

    scaling.factors <- scaled.data %>%
      colMeans()
    scaling.factors <- 1 / scaling.factors

    # multiply each column with corresponding mass
    for (num.factor in seq(1, ncol(scaled.data), 1)) {
      scaled.data[, num.factor] <- scaled.data[, num.factor] * scaling.factors[num.factor]
    }
    scaled.data <- dplyr::bind_cols(G_matrix %>% select(model_type, unit, model_run, date) %>%
      mutate(unit = "normalised"), scaled.data)
    G_matrix <- dplyr::bind_rows(G_matrix, scaled.data)
  } else {
    # check and warn if not scaled to unity
    scaling.factors <- G_matrix %>%
      filter(unit == "normalised") %>%
      select(contains("factor")) %>%
      colMeans()
    all.eqNum <- function(...) as.numeric(sub(".*:", "", all.equal(...)))
    diff <- all.eqNum(mean(scaling.factors), 1)
    if (diff > threshold_unity) {
      cli::cli_warn(c(
        "The normalised {.var G_matrix} is not equal to unity:",
        "i" = "The average contributions of the factors in the normalised
          {.var G_matrix} are not all equal to 1. The difference of {round(diff,2)}
          exceed the value of ({threshold_unity}).",
        "x" = "Consider setting {.var rescale_unity} to `TRUE`."
      ))
    }
  }

  if (!identical(factor_mass, NA)) {
    # Use the mass to calculate concentrations (and if we renormalised, then use that)
    concentration.data <- G_matrix %>%
      filter(unit == "normalised") %>%
      select(contains("factor"))

    # check if the number of columns are equal
    if (ncol(concentration.data) != length(factor_mass)) {
      cli::cli_abort(c(
        "Number of factors must be equal",
        "x" = "The number of factors in {.var G_matrix}
          ({ncol(concentration.data)}) and {.var factor_mass}
          ({length(factor_mass)}) are not equal."
      ))
    }

    # multiply each column with corresponding mass
    for (num.factor in seq(1, ncol(concentration.data), 1)) {
      concentration.data[, num.factor] <- concentration.data[, num.factor] * factor_mass[num.factor]
    }
    ## add additional columns
    concentration.data <- dplyr::bind_cols(G_matrix %>%
                                             select(model_type,
                                                    unit,
                                                    model_run,
                                                    date) %>%
                                             mutate(unit = "concentrations"),
                                           concentration.data)
    ## combine to one dataframe
    G_matrix <- dplyr::bind_rows(G_matrix, concentration.data)
  }

  if(tidy_output) {
    id_variables <- c("model_run", "unit", "date", "model_type")
    # Test for id
    if ("id" %in% names(G_matrix)) {
      id_variables <- c(id_variables, "id")
    }

    if ("run_type" %in% names(G_matrix)) {
      id_variables <- c(id_variables, "run_type")
    }

    G_matrix <- G_matrix %>%
      tidyr::pivot_longer(-dplyr::all_of(id_variables), names_to = "factor") %>%
      group_by(model_run,
               unit,
               date) %>%
      ungroup()


  }
  return(G_matrix)
}
