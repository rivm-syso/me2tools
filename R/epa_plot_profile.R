#' Plot profiles using the EPA-PMF style
#'
#' Function to plot PMF factor profiles using data from PMFR (using
#' \code{\link[pmfr]{read_pmf_factor_profiles}} and
#' \code{\link[pmfr]{tidy_pmf_profiles}}) or from the ME2tools package using
#' \code{\link{me2_read_F}} in the same way as in EPA-PMF.
#'
#' @param F_matrix Tibble from \code{\link[pmfr]{me2_read_F}} containing
#'   the results for the F_matrix. In order to use this function only the
#'   results from a single run should be provided.
#' @param xlabel.angle What angle should the x-axis labels be presented in? If
#'   your labels are long, \code{45} degrees can be useful, which is also the
#'   default.
#' @param xlabel.vjust vertical justification of the xlabel, between \[0,1\].
#'   The default is NA, so that ggplot uses some heuristics to pick the best
#'   value for this parameter. Any other value is processed by changing theme
#'   settings,
#' @param xlabel.hjust horizontal justification of the xlabel, between \[0,1\].
#'   The default is NA, so that ggplot uses some heuristics to pick the best
#'   value for this parameter. Any other value is processed by changing theme
#'   settings,
#' @param xlabel.order The labels containing the species on the x-axis are
#'   plotted based on the factor levels. This parameter can contain the levels
#'   and the species are transformed to factors using these levels. N.B. these
#'   levels can contain code for expressions.
#' @param xlabel.top.bottom Boolean used for splitting the xlabels into a top
#'   and bottom section. The odd indices of the labels will presented at the top
#'   using a secondary axis, whereas the even indices of the labels are
#'   presented at the bottom. This works especially well with plots containing
#'   a lot of species. Instead of using a small font size, the labels are
#'   distributed to the top and bottom. Default setting is \code{FALSE}
#' @param y_min The minimum value for the y-axis, which is used to cut-off the
#'   logarithmic scale. The default setting is \code{-5}, which corresponds to
#'   1xE-5.
#' @param xlab x-axis label.
#' @param ylab Primary y-axis label.
#' @param ylab2 Secondary y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
#' @param x.font.size The size of the xtick labels. Defaults to 10.
#' @param x.n.dodge The number of rows on the x-axis that should be used to 
#'   render the labels. Useful for labels that otherwise overlap, for example
#'   with a large number of species. Also allows for a larger 
#'   \sQuote{x.font.size}.
#' @param bar.color Provide the fill color for the concentration based bars on
#'   the logarithmic scale. Is set to \sQuote{cadetblue3}, as an approximation 
#'   of the color used in EPA-PMF program.
#' @param bar.width The width of the bar, expressed as a value between \[0,1\].
#' @param point.color The fill color for the explained variation, with the
#'   default being \sQuote{firebrick2}, similar to the color used in the EPA-PMF
#'   program.
#' @param point.size The size of the point of the explained variation, 
#'   defaults to 3.
#' @param point.shape The shape of the point of the explained variation. The 
#'   default is 19 (circle). Using 15 will produce a square.
#' @param lollipops Should the dot representing the EV be extended with a line,
#'   representing so it looks more like a lollipop? Defaults to \code{TRUE}.
#' @param errorbar If the F_matrix contains the results of the error estimates,
#'   the error bars for one of those estimates can be plotted in this plot. 
#'   Valid options are \sQuote{none} (default), \sQuote{DISP}, \sQuote{BS} or
#'   \sQuote{BSDISP}.
#' @param errorbar.color The color used for plotting the error bars for
#'   DISP. Default is a list with default colors for the BS, DISP and BSDISP
#'   error bars, based on the colors used in \code{epa_plot_errorsummary}. It
#'   can also be a single color (i.e., \code{errorbar.color = "darkorange3"}. 
#'   If the function cannot match the colors, a default color of 
#'   \sQuote{darkorange3} is applied.
#' @param errorbar.width The width of the upper and lower bars for displaying
#'   the error estimate. Defaults to 0.45 (i.e., the bar.width / 2)
#' @param errorbar.size The line size of the error bars. Defaults to 1.
#' @param errorbar.point.size The size of the point used to display the 
#'   average/median error estimates when \code{cp.run.type} is set as 
#'   \sQuote{base_run}.
#' @param show.plot A logical argument evaluating to TRUE or FALSE indicating
#'   whether the plot should be shown as default.
#' @param expand.mult Vector of multiplicative range expansion factors used on
#'   the x-axis . Defaults to \code{c(0.015,0.005)}, as this seems to work best.
#'   If length is 1, then both left and right x-axis are multiplied with this
#'   factor. If length is 2, the left limit is expanded by \code{expand.mult[1]}
#'   and the right by \code{expand.mult[2]}.
#' @param rm.grid.x Should the vertical grid lines be removed? In some cases
#'   there are a lot of species, causing the vertical grid lines to clutter the
#'   plot. If set to \code{TRUE} these lines are removed. Defaults to
#'   \code{FALSE}.
#' @param perc.x.interval The interval used for the x.ticks on the secondary
#'   axis. The default is 20, leading to intervals of
#'   \dQuote{c(0, 20, 40, 60, 80, 100)}. The value 25 leads to
#'   \dQuote{c(0, 25, 50, 75, 100)}.
#' @param cp.run.type The concentration and percentage of species are typically
#'   plotted using the \sQuote{base_run} values (default). However, by changing
#'   this variable to \sQuote{DISP_avg}, \sQuote{BS_median} or 
#'   \sQuote{BSDISP_avg}, these values (if available) are plotted.
#' @param facet.parse.label Should the labels be parsed using the 
#'   \code{labeller = label_parsed}? If set to \code{TRUE} then \code{"SO[2]"}
#'   will use subscript on the labels shown in the facet.
#' @param ... Other parameters, for example renamed parameters.
#'
#' @return me2tools list containing the ggplot2 with segment and point
#'   geometries and the call.
#'
#' @section Input format of the F-matrix:
#' 
#' The format of the F-matrix as a tibble consists of seven columns. These 
#' seven columns are
#' \itemize{
#'   \item \code{model_type}: a descriptive column with the name of the model 
#'     used to calculate the results (e.g., \dQuote{ME-2} or \dQuote{EPA-PMF})
#'   \item \code{factor_profile}: consists of at least two different 
#'     profiles: \dQuote{concentration_of_species} and 
#'     \dQuote{percentage_of_species_sum}. The concentration profile is in 
#'     concentration units, whereas the percentage is in percentages \[0,100\].
#'   \item \code{model_run}: integer to denote the model run, for example used 
#'     when reading the results of multiple run. When processing only the 
#'     \dQuote{base_run} from RPA-PMF this column can be set to 1.
#'   \item \code{species}: the name of the species
#'   \item \code{run_type}: consists of different types of runs. There should 
#'     be at least one set of data with different \code{factor_profiles} for 
#'     \code{run_type = "base_run"}. Other options for run_type are 
#'     \dQuote{DISP_min}, \dQuote{DISP_avg} and \dQuote{DISP_max} for DISP 
#'     results. For BS results the options are \dQuote{BS_P05}, 
#'     \dQuote{BS_median} and \dQuote{BS_P95}, and for BSDISP results options 
#'     \dQuote{"BSDISP_P05}, \dQuote{BSDISP_avg} and \dQuote{BSDISP_P95} can 
#'     be used. 
#'   \item \code{factor}: the names of the resolved factors.
#'   \item \code{value} : the values, either in concentration units or 
#'     percentages.
#' }
#' 
#' A typical plot of the F-matrix shows the concentration of each species on a 
#' log scale (bars), and the percentage of species sum in a percentage scale 
#' (points). By default, the values associated with both \code{factor_profiles} 
#' \dQuote{concentration_of_species} and \dQuote{percentage_of_species_sum} for 
#' the \code{run_type} \dQuote{base_run} is used.
#' 
#' @section Plotting of errorbars:
#' 
#' The function is capable of plotting errorbars for the results of either DISP, 
#' BS or BSDISP. Only one set of results is supported at this point. To use 
#' this feature, \code{errorbar} should be set to one of the options 
#' \dQuote{DISP}, \dQuote{BS} or \dQuote{BSDISP} (default is \dQuote{none}). 
#' Based on the selection the function automatically checks for and uses the 
#' provided data.
#' 
#' When using \code{cp.run.type = "base_run"}, meaning that the concentration 
#' bars and the percentage dots are plotted using the values associated with 
#' the base_run, plotting the error bars only require the run_type of 
#' \dQuote{DISP_min} and \dQuote{DISP_max} for DISP; the \dQuote{BS_P05} and 
#' \dQuote{BS_P95} for BS and \dQuote{BSDISP_P05} and \dQuote{BSDISP_P95} for 
#' BSDISP, respectively. The average/median of the error estimates is plotted 
#' using \dQuote{DISP_avg}, \dQuote{BS_median}, or \dQuote{BSDISP_avg} and 
#' displayed as a point within the error bars when 
#' \code{cp.run.type = "base_run"}. The selected \code{run_type} should at least
#' be available in concentration units (i.e., 
#' \code{factor_profile = "concentration_of_species"}).
#' 
#' The concentration bars and percentage dots can also be plotted based on the 
#' average/median of the error estimates. In that case it is important that for 
#' each average/median of the error estimates also the percentages of the sum 
#' of the species (i.e., \code{factor_profile = "percentage_of_species_sum"}) 
#' is present in the F-matrix.
#' 
#' @export
#'
#' @import cli
#' @import openair
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#'
epa_plot_profile <- function(F_matrix,
                             xlabel.angle = 45,
                             xlabel.vjust = NA,
                             xlabel.hjust = NA,
                             xlabel.order = NA,
                             xlabel.top.bottom = FALSE,
                             y_min = -5,
                             ylab = "Concentration of species",
                             ylab2 = "EV (%)",
                             xlab = "Species",
                             auto.text = TRUE,
                             x.font.size = 10,
                             x.n.dodge = 1,
                             bar.color = "cadetblue3",
                             bar.width = 0.9,
                             point.color = "firebrick2",
                             point.size = 3,
                             point.shape = 19,
                             lollipops = TRUE,
                             errorbar = "none",
                             errorbar.color = list("BS" = "goldenrod2",
                                                   "BSDISP" = "green4",
                                                   "DISP" = "royalblue2"),
                             errorbar.width = 0.45,
                             errorbar.size = 1,
                             errorbar.point.size = 3,
                             show.plot = TRUE,
                             expand.mult = c(0.015, 0.005),
                             rm.grid.x = FALSE,
                             perc.x.interval = 20,
                             cp.run.type = "base_run",
                             facet.parse.label = FALSE,
                             ...) {

  # The EPA factor profile plot consists of a dual-axis plot containing the
  # log10 concentrations and the percentage of each species in a plot. At first
  # sight the concentrations look like a bar plot, but this is not the case.
  # A bar plot typically starts from 0 (within ggplot this can be changed to
  # any arbitrary number).
  # Instead in this plotting routine we use geometries to plot the bars with a
  # lower limit of 0.00001 (-5)


  ##################################################################
  ##                            Checks                            ##
  ##################################################################

  # check if tidied
  if (!"factor" %in% names(F_matrix)) {
    cli::cli_abort(c(
      "{.var factor} column not detected:",
      "i" = "The F-matrix should contain a {.var factor} column.",
      "x" = "Did you forget to enable {.var tidy_output} when reading the F-matrix?"
    ))
  } else {
    num.factors <- length(unique(F_matrix$factor))
  }

  # check if only one model_run
  if ("model_run" %in% names(F_matrix)) {
    if (length(unique(F_matrix$model_run)) > 1) {
      cli::cli_abort(c(
        "More than 1 run detected:",
        "i" = "The F-matrix contains more than 1 {.var model_run}.",
        "x" = "Did you select a base case run using a filter on {.var model_run}?"
      ))
    }
  }

  # check if errorbar has a valid value
  if (!errorbar %in% c("none", "DISP", "BS", "BSDISP")) {
    cli::cli_abort(c(
      "Unknown value for {.var errorbar}:",
      "i" = "'{errorbar}' is not a valid value for {.var errorbar}.",
      "x" = "Valid values are 'none' (default), 'DISP', 'BS' or 'BSDISP'"
    ))
  }

  # check if cp.run.type has valid value
  if (!cp.run.type %in% c("base_run", "DISP_avg", "BS_median", "BSDISP_avg")) {
    cli::cli_abort(c(
      "Unknown value for {.var cp.run.type}:",
      "i" = "'{cp.run.type}' is not a valid value for {.var cp.run.type}.",
      "x" = "Valid values are 'base_run' (default), 'DISP_avg', 'BS_median' or 'BSDISP_avg'"
    ))
  }

  # Check and adjust labels
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
  if (!identical(ylab2, NA)) {
    if (auto.text) {
      ylab2 <- openair::quickText(ylab2)
    }
  }

  # check if it contains "concentration_of_species" and
  # "percentage_of_species_sum"

  # check if species is factor
  if (!"factor" %in% class(F_matrix$species)) {
    F_matrix <- F_matrix %>%
      mutate(species = factor(paste0("`", species, "`")))
  }

  # create the plot.data
  df <- F_matrix %>%
    dplyr::mutate(value = dplyr::if_else(factor_profile == "concentration_of_species",
      log10(dplyr::if_else(value < 10^y_min, 10^y_min, value)),
      as.double(value)
    )) %>%
    dplyr::mutate_if(is.numeric, list(~ na_if(., Inf))) %>%
    dplyr::mutate_if(is.numeric, list(~ na_if(., -Inf))) %>%
    dplyr::mutate(value = tidyr::replace_na(value, y_min))

  if (nrow(df) == 0) {
    cli::cli_abort(c(
      "No data selected to plot:",
      "i" = "The F-matrix seems to be empty.",
      "x" = "Did you select the correct data for {.var F_matrix}?"
    ))
  }
  
  # check colors
  if (length(bar.color) > 1) {
    if (length(bar.color) != num.factors) {
      cli::cli_abort(c(
        "Not enough bar colors:",
        "i" = "The number of bar colors needs to be equal to 1 or the number of factors.",
        "x" = "Did you provide the correct amount of colors in {.var bar.color}?"
      ))
    }
  }
  
  # errorbar colors
  if (is.list(errorbar.color)) {
    if (exists(errorbar, where=errorbar.color)) {
      errorbar.color <- errorbar.color[[errorbar]]
    } else {
      errorbar.color <- "darkorange3"
    }
  } else {
    if (length(errorbar.color) != 1) {
      errorbar.color <- "darkorange3"
    }
  }

  #################################################################
  ##                    Checks for error data                    ##
  #################################################################

  # prepare DISP data if needed
  # check if we need to do something with the DISP profiles so we prepare the
  # data here
  if (errorbar == "DISP") {
    check_vars <- c("DISP_avg", "DISP_min", "DISP_max")
    for (check_var in check_vars) {
      if (!(check_var %in% unique(df$run_type))) {
        cli::cli_abort(c(
          "DISP results not found:",
          "i" = "The F-matrix does not contain a {.var run_type} column containing {check_var}.",
          "x" = "Did you provide the DISP results with the {.var F_matrix}?"
        ))
      } else {
        # check for "concentration_of_species" in combination with run_type
        test.df <- df %>%
          filter(
            run_type == check_var,
            factor_profile == "concentration_of_species"
          )
        if (nrow(test.df) == 0) {
          cli::cli_abort(c(
            "No concentrations for specific run_type:",
            "i" = "{.var factor_profile} = 'concentration_of_species' does not exists for {.var run_type}={check_var}.",
            "x" = "Did you forget to add the concentration values for DISP?"
          ))
        }
      }
    }
    # check for percentage of species if we need to plot those
    if (cp.run.type == "DISP_avg") {
      test.df <- df %>%
        filter(
          run_type == "DISP_avg",
          factor_profile == "percentage_of_species_sum"
        )
      if (nrow(test.df) == 0) {
        cli::cli_abort(c(
          "No percentages for specific run_type:",
          "i" = "{.var factor_profile} = 'percentage_of_species_sum' does not exists for {.var run_type}='DISP_avg'.",
          "x" = "Did you forget to add the percentage of species sum for DISP?"
        ))
      }
    }
    # set the correct profile for concentrations
    error_run_type <- "DISP_"
    error_ymin <- "DISP_min"
    error_ymax <- "DISP_max"
    error_yavg <- "DISP_avg"
  }

  # prepare BS data if needed
  # check if we need to do something with the DISP profiles so we prepare the
  # data here
  if (errorbar == "BS") {
    check_vars <- c("BS_median", "BS_P05", "BS_P95")
    for (check_var in check_vars) {
      if (!(check_var %in% unique(df$run_type))) {
        cli::cli_abort(c(
          "BS results not found:",
          "i" = "The F-matrix does not contain a {.var run_type} column containing {check_var}.",
          "x" = "Did you provide the BS results with the {.var F_matrix}?"
        ))
      } else {
        # check for "concentration_of_species" in combination with run_type
        test.df <- df %>%
          filter(
            run_type == check_var,
            factor_profile == "concentration_of_species"
          )
        if (nrow(test.df) == 0) {
          cli::cli_abort(c(
            "No concentrations for specific run_type:",
            "i" = "{.var factor_profile} = 'concentration_of_species' does not exists for {.var run_type}={check_var}.",
            "x" = "Did you forget to add the concentration values for BS?"
          ))
        }
      }
    }
    # check for percentage of species if we need to plot those
    if (cp.run.type == "BS_median") {
      test.df <- df %>%
        filter(
          run_type == "BS_median",
          factor_profile == "percentage_of_species_sum"
        )
      if (nrow(test.df) == 0) {
        cli::cli_abort(c(
          "No percentages for specific run_type:",
          "i" = "{.var factor_profile} = 'percentage_of_species_sum' does not exists for {.var run_type}='BS_median'.",
          "x" = "Did you forget to add the percentage of species sum for DISP?"
        ))
      }
    }
    # set the correct profile for concentrations
    error_run_type <- "BS_"
    error_ymin <- "BS_P05"
    error_ymax <- "BS_P95"
    error_yavg <- "BS_median"
  }

  # prepare BSDISP data if needed
  # check if we need to do something with the DISP profiles so we prepare the
  # data here
  if (errorbar == "BSDISP") {
    check_vars <- c("BSDISP_avg", "BSDISP_P05", "BSDISP_P95")
    for (check_var in check_vars) {
      if (!(check_var %in% unique(df$run_type))) {
        cli::cli_abort(c(
          "BSDISP results not found:",
          "i" = "The F-matrix does not contain a {.var run_type} column containing {check_var}.",
          "x" = "Did you provide the BSDISP results with the {.var F_matrix}?"
        ))
      } else {
        # check for "concentration_of_species" in combination with run_type
        test.df <- df %>%
          filter(
            run_type == check_var,
            factor_profile == "concentration_of_species"
          )
        if (nrow(test.df) == 0) {
          cli::cli_abort(c(
            "No concentrations for specific run_type:",
            "i" = "{.var factor_profile} = 'concentration_of_species' does not exists for {.var run_type}={check_var}.",
            "x" = "Did you forget to add the concentration values for BSDISP?"
          ))
        }
      }
    }
    # check for percentage of species if we need to plot those
    if (cp.run.type == "BSDISP_avg") {
      test.df <- df %>%
        filter(
          run_type == "BSDISP_avg",
          factor_profile == "percentage_of_species_sum"
        )
      if (nrow(test.df) == 0) {
        cli::cli_abort(c(
          "No percentages for specific run_type:",
          "i" = "{.var factor_profile} = 'percentage_of_species_sum' does not exists for {.var run_type}='BSDISP_avg'.",
          "x" = "Did you forget to add the percentage of species sum for DISP?"
        ))
      }
    }
    # set the correct profile for concentrations
    error_run_type <- "BSDISP_"
    error_ymin <- "BSDISP_P05"
    error_ymax <- "BSDISP_P95"
    error_yavg <- "BSDISP_avg"
  }

  #################################################################
  ##             Plot preparations (axis breaks etc)             ##
  #################################################################

  # order the x_labels
  if (!identical(xlabel.order, NA)) {
    df$species <- factor(as.character(df$species), levels = xlabel.order)
  }

  # create numeric column based on species. We use this to plot the data.
  df$numeric_x <- as.numeric(df$species)

  # Calculate the upper axis value
  max_conc_log <- df %>%
    filter(
      factor_profile == "concentration_of_species",
      stringr::str_detect(run_type, cp.run.type)
    ) %>%
    select(value) %>%
    max()

  ## calculate the optimum axis
  breaks_log <- seq(y_min, ceiling(max_conc_log), 2)

  # if the breaks are to small, then increase them to the minimum
  # special situation
  if ((length(breaks_log) == 4) && (max(breaks_log = 1))) {
    breaks_log <- c(breaks_log, max(breaks_log) + 2)
  }

  # bare minimum is 4
  if ((length(breaks_log) < 4)) {
    breaks_log <- c(breaks_log, max(breaks_log) + 2)
  }

  labels_log <- paste0("10^", breaks_log)

  breaks_percentage <- seq(y_min, max(breaks_log), 1) * (length(seq(y_min, max(breaks_log), 1)) - 1) + y_min

  # calculate the gap
  gap_percentage <- (max(breaks_percentage) - min(breaks_percentage)) / (100 / perc.x.interval)
  sec.breaks <- c(min(breaks_percentage), round((rep(1, (100 / perc.x.interval)) * min(breaks_percentage)) + (seq(1, (100 / perc.x.interval), 1) * gap_percentage), 0))

  # scale_items <- ceiling(length(breaks_percentage) / 2)
  # sec.breaks <- c(seq(min(breaks_percentage), max(breaks_percentage), round(length(seq(min(breaks_percentage), max(breaks_percentage), 1)) / scale_items, 0)), max(breaks_percentage))
  sec.labels <- round(seq(0, 100, (100 / (length(sec.breaks) - 1))), 1)

  # scale the percentages to the secondary axis
  df <- df %>%
    dplyr::mutate(value = if_else(factor_profile == "percentage_of_species_sum", (value / 100) * (length(seq(y_min, max(breaks_log), 1)) - 1) + y_min, as.double(value)))

  
  ##################################################################
  ##                          Set colors                          ##
  ##################################################################
  plot.df <- df %>%
    filter(
      factor_profile == "concentration_of_species",
      stringr::str_detect(run_type, cp.run.type)
    )
  
  if (class(F_matrix$factor) == "factor") {
    myColors <- tibble("factor" = levels(F_matrix$factor),
                       "color" = openair::openColours(bar.color, num.factors))
  } else {
    myColors <- tibble("factor" = unique(F_matrix$factor),
                       "color" = openair::openColours(bar.color, num.factors))
  }
  
  myColors <- left_join(plot.df %>% select(factor),
            myColors,
            by = "factor")
  
  myColors <- myColors[["color"]]
  
  #################################################################
  ##                         Create plot                         ##
  #################################################################

  plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = plot.df,
      ggplot2::aes(
        xmin = numeric_x - (bar.width / 2),
        xmax = numeric_x + (bar.width / 2),
        ymin = y_min, ymax = value
      ),
      col = myColors,
      fill = myColors
    ) +
    ggplot2::geom_point(ggplot2::aes(x = numeric_x, y = value),
      data = df %>%
        filter(
          factor_profile == "percentage_of_species_sum",
          stringr::str_detect(run_type, cp.run.type)
        ),
      size = point.size,
      shape = point.shape,
      col = point.color
    )

  # plot lollipop lines
  if (lollipops) {
    plot <- plot +
      ggplot2::geom_segment(
        data = df %>%
          filter(
            factor_profile == "percentage_of_species_sum",
            stringr::str_detect(run_type, cp.run.type)
          ),
        ggplot2::aes(
          x = numeric_x,
          xend = numeric_x,
          y = y_min,
          yend = value
        ),
        color = point.color,
        alpha = 0.5
      )
  }

  # plot errorbars
  if (errorbar != "none") {
    disp.df <- df %>%
      filter(
        factor_profile == "concentration_of_species",
        stringr::str_detect(run_type, error_run_type)
      ) %>%
      select(numeric_x, factor, run_type, value) %>%
      pivot_wider(names_from = run_type, values_from = value)

    plot <- plot +
      geom_errorbar(
        data = disp.df, aes(
          x = numeric_x,
          ymin = !!sym(error_ymin),
          ymax = !!sym(error_ymax)
        ),
        width = errorbar.width,
        color = errorbar.color,
        size = errorbar.size,
        position = position_dodge(.9)
      )

    if (cp.run.type == "base_run") {
      # add point for the average/median
      plot <- plot +
        ggplot2::geom_point(ggplot2::aes(x = numeric_x, y = !!sym(error_yavg)),
          data = disp.df,
          size = errorbar.point.size,
          shape = 18,
          col = errorbar.color,
          position = position_dodge(.9)
        )
    }
  }
  
  if (facet.parse.label) {
    plot <- plot +
      ggplot2::facet_grid(factor ~ ., labeller = label_parsed)  
    
  } else {
    plot <- plot +
      ggplot2::facet_grid(factor ~ .)  
  }

  plot <- plot +
    ggplot2::ylab(ylab) +
    ggplot2::scale_y_continuous(
      limits = c(y_min, max(breaks_log)),
      breaks = breaks_log,
      labels = parse(text = labels_log),
      sec.axis = ggplot2::sec_axis(
        trans = ~ . * (length(seq(y_min, max(breaks_log), 1)) - 1) + y_min,
        name = ylab2,
        breaks = sec.breaks,
        labels = sec.labels
      ),
      expand = expansion(mult = c(0, 0.05))
    ) +
    ggplot2::annotation_logticks(sides = "l") +
    ggplot2::theme_bw()
  
  # get the labels
  x_labels <- levels(df$species)
  
  # split into even and odd
  x_labels_top <- x_labels
  x_labels_bottom <- x_labels
  # define quick functions
  odd <- function(x) x%%2 != 0
  even <- function(x) x%%2 == 0
  # split the labels
  x_labels_top[odd(seq(1, length(x_labels_top)))] <- NA
  x_labels_bottom[even(seq(1, length(x_labels_bottom)))] <- NA
  
  if ((identical(xlabel.vjust, NA)) & (identical(xlabel.hjust, NA))) {
    # set parameters using guide_axis
    if (xlabel.top.bottom) {
      plot <- plot +
        ggplot2::scale_x_continuous(
          name = xlab,
          labels = parse(text = x_labels_top),
          breaks = sort(unique(df$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = guide_axis(angle = xlabel.angle, n.dodge = x.n.dodge),
          sec.axis = ggplot2::sec_axis(
            trans = ~ .,
            breaks = sort(unique(df$numeric_x)),
            labels =  parse(text = x_labels_bottom),
            guide = ggplot2::guide_axis(angle = xlabel.angle, n.dodge = x.n.dodge)
          )
        ) + 
        theme(axis.text.x.top = ggplot2::element_text(
          hjust = 0,
          vjust = 0.5
        ))
    } else {
      plot <- plot +
        ggplot2::scale_x_continuous(
          name = xlab,
          labels = parse(text = x_labels),
          breaks = sort(unique(df$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = ggplot2::guide_axis(angle = xlabel.angle, n.dodge = x.n.dodge)
        )
    }
    plot <- plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          size = x.font.size
        ),
        axis.ticks.y.right = ggplot2::element_line(color = point.color),
        axis.line.y.right = ggplot2::element_line(color = point.color),
        axis.text.y.right = ggplot2::element_text(color = point.color),
        axis.title.y.right = ggplot2::element_text(color = point.color),
        panel.grid.major.x = element_blank()
      )  
  } else {
    # set parameters using theme
    if (xlabel.top.bottom) {
      plot <- plot +
        ggplot2::scale_x_continuous(
          name = xlab,
          labels = parse(text = x_labels_top),
          breaks = sort(unique(df$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = guide_axis(n.dodge = x.n.dodge),
          sec.axis = ggplot2::sec_axis(
            trans = ~ .,
            breaks = sort(unique(df$numeric_x)),
            labels =  parse(text = x_labels_bottom),
            guide = ggplot2::guide_axis(n.dodge = x.n.dodge)
          )
        ) + 
        theme(axis.text.x.top = ggplot2::element_text(
          hjust = 0,
          vjust = 0.5
        ))
    } else {
      plot <- plot +
        ggplot2::scale_x_continuous(
          name = xlab,
          labels = parse(text = x_labels),
          breaks = sort(unique(df$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = ggplot2::guide_axis(n.dodge = x.n.dodge)
        )
    }
    plot <- plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = xlabel.angle,
          size = x.font.size,
          hjust = xlabel.hjust,
          vjust = xlabel.vjust
        ),
        axis.ticks.y.right = ggplot2::element_line(color = point.color),
        axis.line.y.right = ggplot2::element_line(color = point.color),
        axis.text.y.right = ggplot2::element_text(color = point.color),
        axis.title.y.right = ggplot2::element_text(color = point.color),
        panel.grid.major.x = element_blank()
      )
  }
  
  # remove vertical grid lines?
  if (rm.grid.x) {
    plot <- plot +
      ggplot2::theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  }

  # run garbage collector
  gc()
  if (!identical(dev.list(), NULL)) {
    dev.off()
  }

  ##################################################################
  ##                        Prepare output                        ##
  ##################################################################

  # print(metcor.plot)
  output <- list(
    "plot" = plot,
    call = match.call()
  )
  class(output) <- "me2tools"
  if (show.plot) {
    plot(output$plot)
  }
  invisible(output)
}
