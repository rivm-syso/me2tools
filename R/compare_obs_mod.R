#' Compare the observed total mass vs the modelled contributions in the G-matrix
#'
#' Based on the observed total mass and the modelled contributions in the
#' G-matrix this function performs a multilinear regression on the observations
#' and the individual normalised factor contributions. The function also
#' performs a regression against the observations and the sum of the factors
#'
#' @param measurements Tibble containing the date and the observations against
#'   which a regression should be performed. Typical usage would be to provide
#'   the total mass (i.e., PM2.5 or PM10).
#' @param obs The name of the column in the "measurement" tibble containing the
#'   observations for the regression.
#' @param G_matrix The preliminary output from functions reading G values in a
#'   non-tidied format.
#' @param mod Vector containing the names of the columns in the G matrix which
#'   are considered the factors and should be used in the regression.
#' @param factor_mass This is a vector with the same length as the number of
#'   factor profiles. It contains the total mass in concentration units which is
#'   used to transform the G matrix from unity to concentration units.
#' @param robust Boolean used to switch between normal (OLS) or robust (RLS)
#'   regression.
#' @param show.plot A logical argument evaluating to TRUE or FALSE indicating
#'   whether the plot should be shown as default.
#' @param line.color The color of the regression line. Defaults to
#'   "dodgerblue4".
#' @param line.size The thickness of the regression line. Defaults to 1.
#' @param point.color The color of the points. Defaults to "red".
#' @param point.shape The shape of the points. Defaults to 16.
#' @param point.size The size of the points. Defaults to 2.
#' @param point.alpha The alpha transparency for the points. Defaults to 0.25.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
#'
#' @returns me2tools list containing the ggplot2 for the regression of the
#'  observations and the sum of the factors in concentration units (LR). The
#'  results of both the MLR and the LR models are also given, as well as the
#'  data used to calculate these models. Additionally information is provided
#'  regarding the use of the concentration units in G and if robust regression
#'  has been applied.
#'
#'  @section Using factor_mass for the sum of factors regression:
#'  To perform the regression against the sum of factors, the factors in the
#'  G-matrix also have to be available in the concentration units. By default,
#'  only the "normalised" data is available, and a \code{factor_mass} is used to
#'  transform this data into concentration units.
#'
#'  If the provided G-matrix only contains normalised data, a vector of the same
#'  length of the number of factors can be provided in this function to
#'  calculate the concentration units.
#'
#'  If this vector, called \code{factor_mass}, is not provided and no
#'  concentration units are present in G, then the coefficients of the MLR are
#'  used as the \code{factor_mass}.
#'
#'  If concentration units are present in G, than this will take precedence,
#'  unless the \code{factor_mass} is defined. In the latter case the existing
#'  concentration data is overwritten with the new concentrations based on the
#'  provided \code{factor_mass}.
#'
#'  In the output of the function a variable is used to denote the application
#'  of \code{factor_mass} and the possible concentrations in G matrix in the
#'  regression.
#'
#' @import cli
#' @importFrom MASS rlm
#' @importFrom stats lm
#' @import dplyr
#' @import ggplot2
#' @importFrom openair quickText
#'
compare_obs_mod <- function(measurements,
                            obs = "pm10",
                            G_matrix,
                            mod = c(
                              "factor_01",
                              "factor_02"
                            ),
                            factor_mass = NA,
                            robust = FALSE,
                            show.plot = TRUE,
                            line.color = "dodgerblue4",
                            line.size = 1,
                            point.color = "red",
                            point.shape = 16,
                            point.size = 2,
                            point.alpha = 0.25,
                            ylab = "Modelled",
                            xlab = "Observed",
                            auto.text = TRUE) {
  if (!"normalised" %in% unique(G_matrix$unit)) {
    cli::cli_abort(c(
      "{.var unit} == 'normalised' not detected:",
      "i" = "The {.var G_matrix} should provide rows with {.var unit} = 'normalised'.",
      "x" = "Did you provide the correct {.var G_matrix} to be used?"
    ))
  }

  # Check if measurements has the same length as G_matrix
  if (nrow(measurements) != nrow(G_matrix %>% filter(unit == "normalised"))) {
    # error size
    cli::cli_abort(c(
      "Different dimensions:",
      "i" = "The G-matrix and measurements have different number of rows.",
      "x" = "Did you provide the correct data?"
    ))
  }

  # check dates
  single_G <- G_matrix %>% filter(unit == "normalised")
  date_diff <- measurements$date - single_G$date
  if (as.numeric(min(date_diff)) != 0 | as.numeric(max(date_diff)) != 0) {
    # difference in dates
    cli::cli_abort(c(
      "Dates do not match:",
      "i" = "The G-matrix and measurements have different dates.",
      "x" = "Did you provide the correct data?"
    ))
  }
  # check if the mod variables exists
  for (factor in mod) {
    if (!factor %in% names(G_matrix)) {
      cli::cli_abort(c(
        "{.var factor} column not detected:",
        "i" = "The G-matrix should contain a {.var factor} column.",
        "x" = "Did you provide the correct {.var mod} to be used?"
      ))
    }
  }

  if (!obs %in% names(measurements)) {
    cli::cli_abort(c(
      "{.var obs} column not detected:",
      "i" = "The measurements should contain a {obs} column.",
      "x" = "Did you provide the correct {.var obs} to be used?"
    ))
  }

  # check mass factors
  if (length(factor_mass) == 1) {
    if (!is.na(factor_mass)) {
      cli::cli_abort(c(
        "Length should be equal to the number of factors:",
        "i" = "{.var factor_mass} has length {length(factor_mass)}.",
        "x" = "Length should be equal to {length(mod)}"
      ))
    }
  } else {
    if (length(factor_mass) != length(mod)) {
      cli::cli_abort(c(
        "Length should be equal to the number of factors:",
        "i" = "{.var factor_mass} has length {length(factor_mass)}.",
        "x" = "Length should be equal to {length(mod)}"
      ))
    }
  }

  # join measurements with scaled contributions
  mlr_data <- left_join(measurements,
    G_matrix %>%
      filter(unit == "normalised") %>%
      select(all_of(c("date", mod))),
    by = "date"
  )

  # set labels using quickText
  if (!identical(xlab, NA)) {
    if (auto.text) {
      xlab <- openair::quickText(xlab)
    }
  }
  if (!identical(ylab, NA)) {
    if (auto.text) {
      ylab <- openair::quickText(ylab)
    }
  }

  ##############################################################################
  # see if we can create a model using variables
  fm <- as.formula(paste(obs, "~ 0 +", paste(mod, collapse = " + ")))
  if (robust) {
    modelMLR <- MASS::rlm(fm, data = mlr_data, maxit = 200)
  } else {
    modelMLR <- stats::lm(fm, data = mlr_data)
  }

  ##############################################################################
  if (identical(factor_mass, NA)) {
    if (!"concentrations" %in% unique(G_matrix$unit)) {
      # there are no concentrations in the G_matrix, we now use the coefficients
      # of the MLR as the factor mass
      factor_mass <- as.numeric(modelMLR$coefficients)
      used_mass <- "Calculated MLR coefficients have been used to calculate G contributions in concentrations units."
    } else {
      used_mass <- "Provided concentrations in G have been used as concentration units."
    }
  } else {
    used_mass <- "Provided factor_mass has been used to calculate G contributions in concentrations units."
  }

  ##############################################################################
  # apply factor_mass
  if (!identical(factor_mass, NA)) {
    # Use the mass to calculate concentrations (and if we renormalised, then use that)
    concentration.data <- G_matrix %>%
      filter(unit == "normalised") %>%
      select(all_of(mod))

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
    concentration.data <- dplyr::bind_cols(
      G_matrix %>%
        select(
          model_type,
          unit,
          model_run,
          date
        ) %>%
        mutate(unit = "concentrations"),
      concentration.data
    )
    ## combine to one dataframe
    G_matrix <- dplyr::bind_rows(
      G_matrix %>%
        filter(unit == "normalised"),
      concentration.data
    )
  }

  # create the data for regression
  regr_data <- left_join(measurements,
    G_matrix %>%
      filter(unit == "concentrations") %>%
      mutate(pred = rowSums(G_matrix %>%
        filter(unit == "concentrations") %>%
        select(all_of(mod)))) %>%
      select(all_of(c("date", "pred"))),
    by = "date"
  )

  ##############################################################################
  # see if we can create a model using variables
  fm <- as.formula(paste(obs, "~ 0 + pred"))
  if (robust) {
    modelLM <- MASS::rlm(fm, data = regr_data, maxit = 200)
  } else {
    modelLM <- stats::lm(fm, data = regr_data)
  }

  # create plot
  scatter <- ggplot2::ggplot(regr_data,
                             ggplot2::aes(x = !!sym(obs), y = pred)) +
    geom_point(shape = point.shape,
               colour = point.color,
               alpha = point.alpha,
               size = point.size)
  # specific options for robust or regular regression
  if (robust) {
    scatter <- scatter +
      ggplot2::geom_smooth(
        method = MASS::rlm,
        se = FALSE,
        color = line.color,
        size = line.size,
        formula = y ~ x + 0
      )
  } else {
    scatter <- scatter +
      ggplot2::geom_smooth(
        method = "lm",
        se = FALSE,
        color = line.color,
        size = line.size,
        formula = y ~ x + 0
      )
  }
  xmax <- ggplot_build(scatter)$layout$panel_params[[1]]$x_range[2]
  ymax <- ggplot_build(scatter)$layout$panel_params[[1]]$y_range[2]

  max_axis <- max(c(
    ceiling(ggplot_build(scatter)$layout$panel_scales_x[[1]]$range$range / 10) * 10,
    ceiling(ggplot_build(scatter)$layout$panel_scales_y[[1]]$range$range / 10) * 10
  ))

  scatter <- scatter +
    ggplot2::geom_abline(intercept = 0,
                         slope = 1,
                         color = "black",
                         linetype = "solid",
                         size = 0.25) +
    ggplot2::geom_abline(intercept = 0,
                         slope = 0.5,
                         color = "black",
                         linetype = "dashed",
                         size = 0.25) +
    ggplot2::geom_abline(intercept = 0,
                         slope = 2,
                         color = "black",
                         linetype = "dashed",
                         size = 0.25) +
    ggplot2::theme_bw() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = c(0, max_axis),
                         ylim = c(0, max_axis),
                         expand = FALSE,
                         clip = "on") +
    xlab(xlab) +
    ylab(ylab)



  #################
  # output
  #################
  # print(metcor.plot)
  output <- list(
    "plot" = scatter,
    "mass_used" = used_mass,
    "model_MLR" = modelMLR,
    "data_MLR" = mlr_data,
    "model_LR" = modelLM,
    "data_LR" = regr_data,
    "robust_regression" = robust,
    call = match.call()
  )
  class(output) <- "me2tools"
  if (show.plot) {
    plot(output$plot)
  }
  invisible(output)
}
