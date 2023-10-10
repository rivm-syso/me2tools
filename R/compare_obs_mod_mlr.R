#' Compare the observed total mass vs the modelled contributions in the G-matrix
#'
#' Based on the observed total mass and the modelled contributions in the
#' G-matrix this function performs a multilinear regression (MLR) on the 
#' observations and the individual normalised factor contributions. From this 
#' regression the coefficient of each factor shows the contribution of these 
#' factors to the total mass. If the variable \code{factor_mass} is not 
#' provided, then the results of the MLR is used to transform the normalized 
#' factor contributions towards concentration contributions. The sum of the 
#' contributions of factors in concentration units, either calculated in the 
#' function or provided in the function call, are then used in a regression 
#' against the measured mass. For usage of multilinear regression (MLR) see 
#' function \code{\link{compare_obs_mod_mlr}}
#' 
#' @param measurements Tibble containing the date and the observations against
#'   which a regression should be performed. Typical usage would be to provide
#'   a tibble containing the date and the total mass (i.e., PM2.5 or PM10).
#' @param obs The name of the column in the "measurement" tibble containing the
#'   observations for the regression.
#' @param G_matrix The preliminary output from functions reading G values in a
#'   non-tidied format. In order to perform MLR this matrix needs to have
#'   "normalized" factor contributions. If this matrix also contains the factor
#'   contributions in "concentration" units, then these results are used to
#'   calculate the sum of factors. This sum will then be used in the final 
#'   regression agains the measured mass. Note that this will note use the output
#'   from the MLR nor the \code{factor_mass} to calculate factor contributions
#'   in concentration units based on these results (i.e., provided factor 
#'   contributions always take precedence).
#' @param mod Character vector containing the names of the columns in the G 
#'   matrix which are considered the factors and should be used in the regression.
#' @param factor_mass This is a vector with the same length as the number of
#'   factor profiles. It contains the total mass in concentration units which is
#'   used to transform the G matrix from unity to concentration units. Providing
#'   this factor will take precedence over the MLR results. If the MLR results
#'   should be used, then this factor should be set to \code{NA}, which is also
#'   the default
#' @param robust Boolean used to switch between normal (OLS) or robust (RLS)
#'   regression for both the MLR as the LR.
#' @param maxit The limit on the number of IWLS iterations used when robust
#'   regression is performed. Defaults to 200. Can be increased in case of error
#'   messages related to failed convergence in the number of steps used.
#' @param intercept Should an intercept be used when doing the LR regression
#'   (i.e., the sum of factors vs measured mass)? Default is FALSE, forcing 
#'   the regression through the origin. Note that the MLR is forced through the
#'   origin by default.
#' @param regress.details Show the equation and R2 (non-robust only) on the 
#'   LR plot. Default usage is \code{TRUE} or \code{FALSE}. 
#' @param mod.line This option will add three lines to assist with model
#'   evaluations. These three lines are the 1:1 line (solid) and the 1:0.5 and
#'   1:2 lines (dashed). Adding these lines will help to evaluate how close the
#'   points are to the 1:1 relation. Points within the dashed lines will be
#'   within a factor of two of the regression line.
#' @param CI Should the confidence levels (CI) be shown? Default of this 
#'   parameter is set to \dQuote{FALSE}.
#' @param CI.type How should the CI be displayed? There are two options: 
#'   \dQuote{area} and \dQuote{line}. With \dQuote{area} a shaded area is 
#'   created, whereas \dQuote{line} just plots the upper and lower limits of the
#'   interval.
#' @param CI.level = Confidence level, between \[0,1\], defaults to 
#'   \dQuote{0.95} corresponding to 95%.
#' @param CI.color The color of the CI area or line, depending on the 
#'   \code{CI.type}. Defaults to "dodgerblue4".
#' @param CI.linetype The type of the lines. Is ignored when 
#'   \code{CI.type = "area"}
#' @param CI.linewidth The width of the lines. Is ignored when 
#'   \code{CI.type = "area"}
#' @param CI.alpha The alpha transparency for the CI area. Is ignored when
#'   \code{PI.type = "area"}
#' @param PI = Should the prediction levels (PI) be shown? Default of this 
#'   parameter is set to \dQuote{FALSE}.
#' @param PI.type How should the PI be displayed? See \code{CI.type} for more
#'   info.
#' @param PI.level = Confidence level, between \[0,1\], defaults to 
#'   \dQuote{0.95} corresponding to 95%.
#' @param PI.color = The color of the PI area or line, depending on the 
#'   \code{PI.type}. Defaults to "dodgerblue4".
#' @param PI.linetype The type of the lines. Is ignored when 
#'   \code{PI.type = "area"}
#' @param PI.linewidth The width of the lines. Is ignored when 
#'   \code{PI.type = "area"}
#' @param PI.alpha =  The alpha transparency for the PI area. Is ignored when
#'   \code{PI.type = "area"} 
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
#' @param asp.ratio.1 Output the plot with an aspect ratio of 1. For this to
#'   work the minimum and maximum of x and y axis are used.
#'
#' @returns me2tools list containing the ggplot2 for the regression of the
#'  observations and the sum of the factors in concentration units (LR). The
#'  results of both the MLR and the LR models are also given, as well as the
#'  data used to calculate these models. Additional information is provided
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
#'  If concentration units are present in G, then this will take precedence,
#'  unless the \code{factor_mass} is defined. In the latter case, the existing
#'  concentration data is overwritten with the new concentrations based on the
#'  provided \code{factor_mass}.
#'
#'  In the output of the function, a variable is used to denote the application
#'  of \code{factor_mass} and the possible concentrations in G matrix in the
#'  regression.
#'
#' @seealso \code{\link{compare_obs_mod}}
#'
#' @import cli
#' @importFrom MASS rlm
#' @importFrom stats lm
#' @import dplyr
#' @import ggplot2
#' @importFrom openair quickText
#' @importFrom stats predict
#' @importFrom ggpmisc stat_poly_eq
#'
compare_obs_mod_mlr <- function(measurements,
                            obs = "pm10",
                            G_matrix,
                            mod = c(
                              "factor_01",
                              "factor_02"
                            ),
                            factor_mass = NA,
                            robust = FALSE,
                            maxit = 200,
                            intercept = FALSE,
                            mod.line = FALSE,
                            regress.details = TRUE,
                            CI = TRUE,
                            CI.type = "area",
                            CI.level = 0.95,
                            CI.color = "dodgerblue4",
                            CI.alpha = 0.4,
                            CI.linetype = "dotted",
                            CI.linewidth = 1,
                            PI = TRUE,
                            PI.type = "area",
                            PI.level = 0.95,
                            PI.color = "dodgerblue4",
                            PI.linetype = "dashed",
                            PI.linewidth = 1,
                            PI.alpha = 0.4,
                            show.plot = TRUE,
                            line.color = "dodgerblue4",
                            line.size = 1,
                            point.color = "red",
                            point.shape = 16,
                            point.size = 2,
                            point.alpha = 0.8,
                            ylab = "Modelled",
                            xlab = "Observed",
                            auto.text = TRUE,
                            asp.ratio.1 = FALSE  ) {
  
  if (!"normalised" %in% unique(G_matrix$unit)) {
    cli::cli_abort(c(
      "{.var unit} == 'normalised' not detected:",
      "i" = "The {.var G_matrix} should provide rows with {.var unit} = 'normalised'.",
      "x" = "Did you provide the correct {.var G_matrix} to be used?"
    ))
  }

  fixed_G = TRUE
  # Check if measurements has the same length as G_matrix
  if (nrow(measurements) != nrow(G_matrix %>% filter(unit == "normalised"))) {
    if (!mod[[1]] %in% names(G_matrix)) {
      # check if data happens to be in the longer format
      if ("factor" %in% names(G_matrix)) {
        id_cols <- c("model_type",
                     "unit",
                     "model_run",
                     "run_type",
                     "date")
        G_matrix <- G_matrix %>% 
          pivot_wider(id_cols = all_of(id_cols),
                      names_from = "factor",
                      values_from = "value")
        fixed_G <- TRUE
      } else {    
        # error size
        cli::cli_abort(c(
          "Different dimensions:",
          "i" = "The G-matrix and measurements have different number of rows.",
          "x" = "Did you provide the correct data?"
        ))
      }
      if (!fixed_G) {
        # error size
        cli::cli_abort(c(
          "Different dimensions:",
          "i" = "The G-matrix and measurements have different number of rows.",
          "x" = "Did you provide the correct data?"
        ))
      }
    }
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
      "i" = "The {.varmeasurements} should contain a {obs} column.",
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

 
  ##############################################################################
  # see if we can create a model using variables
  fm_mlr <- stats::as.formula(paste(obs, "~ 0 +", paste(mod, collapse = " + ")))
  
  if (robust) {
    modelMLR <- MASS::rlm(fm_mlr, data = mlr_data, maxit = 200)
  } else {
    modelMLR <- stats::lm(fm_mlr, data = mlr_data)
  }

  ##############################################################################
  if (identical(factor_mass, NA)) {
    if (!"concentration" %in% unique(G_matrix$unit)) {
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
        mutate(unit = "concentration"),
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
                           filter(unit == "concentration") %>%
                           mutate(pred = rowSums(G_matrix %>%
                                                   filter(unit == "concentration") %>%
                                                   select(all_of(mod)))) %>%
                           select(all_of(c("date", "pred"))),
                         by = "date"
  ) %>% 
    rename(x := !!sym(obs),
           y = pred)
  
  # we got an x and an y variable now, with pred being the sum of the factors
  # Pass these results to the compare_obs_mod function.
 
  model_output <- compare_obs_mod(data = regr_data,
                                  x = "x",
                                  y = "y",
                                  robust = robust,
                                  maxit = maxit,
                                  intercept = intercept,
                                  mod.line = mod.line,
                                  regress.details = regress.details,
                                  CI = CI,
                                  CI.type = CI.type,
                                  CI.level = CI.level,
                                  CI.color = CI.color,
                                  CI.alpha = CI.alpha,
                                  CI.linetype = CI.linetype,
                                  CI.linewidth = CI.linewidth,
                                  PI = PI,
                                  PI.type = PI.type,
                                  PI.level = PI.level,
                                  PI.color = PI.color,
                                  PI.linetype = PI.linetype,
                                  PI.linewidth = PI.linewidth,
                                  PI.alpha = PI.alpha,
                                  show.plot = FALSE,
                                  line.color = line.color,
                                  line.size = line.size,
                                  point.color = point.color,
                                  point.shape = point.shape,
                                  point.size = point.size,
                                  point.alpha = point.alpha,
                                  ylab = ylab,
                                  xlab = xlab,
                                  auto.text = auto.text,
                                  asp.ratio.1 = asp.ratio.1  )
 
  #################
  # output
  #################
  output <- list(
    "plot" = model_output$plot,
    "mass_used" = used_mass,
    "model_MLR" = modelMLR,
    "data_MLR" = mlr_data,
    "formula_LR" = model_output$formula_LR,
    "model_LR" = model_output$model_LR,
    "data_LR" = model_output$data_LR,
    "robust_regression" = model_output$robust_regression,
    "regression.output" = model_output$regression.output,
    call = match.call()
  )
  class(output) <- "me2tools"
  if (show.plot) {
    plot(output$plot)
  }
  invisible(output)
}
