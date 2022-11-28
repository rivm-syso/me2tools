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
#' @param xlabel.hjust horizontal justification of the xlabel, between \[0,1\].
#' @param xlabel.order The labels containing the species on the x-axis are
#'   plotted based on the factor levels. This parameter can contain the levels
#'   and the species are transformed to factors using these levels. N.B. these
#'   levels can contain code for expressions.
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
#' @param bar.color Provide the fill color for the concentration based bars on
#'   the logarithmic scale. Is set to \sQuote{steelblue}, as an approximation of
#'   the color used in EPA-PMF program.
#' @param bar.width The width of the bar, expressed as a value between 0-1.
#' @param point.color The fill color for the explained variation, with the
#'   default being \sQuote{red}, similar to the color used in the EPA-PMF
#'   program.
#' @param point.size The size of the point, typically the bar.width / 2
#' @param lollipops Should the dot representing the EV be extended with a line,
#'   representing so it looks more like a lollipop? Defaults to \code{TRUE}.
#' @param disp.errorbar If the F_matrix contains the DISP min and max values,
#'   the error bars can be plotted in this plot. Note that with this option the
#'   species concentration is replaced by the DISP average.
#' @param disp.errorbar.color The color used for plotting the error bars for
#'   DISP.
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
#'
#' @return me2tools list containing the ggplot2 with segment and point
#'   geometries and the call.
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
                             xlabel.vjust = 0.5,
                             xlabel.hjust = 0.5,
                             xlabel.order = NA,
                             y_min = -5,
                             ylab = "Concentration of species",
                             ylab2 = "EV (%)",
                             xlab = "Species",
                             auto.text = TRUE,
                             x.font.size = 10,
                             bar.color = "cadetblue3",
                             bar.width = 0.9,
                             point.color = "red",
                             point.size = 3,
                             lollipops = TRUE,
                             disp.errorbar = FALSE,
                             disp.errorbar.color = "orange",
                             show.plot = TRUE,
                             expand.mult = c(0.015,0.005),
                             rm.grid.x = FALSE,
                             ...) {

  # The EPA factor profile plot consists of a dual-axis plot containing the
  # log10 concentrations and the percentage of each species in a plot. At first
  # sight the concentrations look like a bar plot, but this is not the case.
  # A bar plot typically starts from 0 (within ggplot this can be changed to
  # any arbitrary number).
  # Instead in this plotting routine we use geometries to plot the bars with a
  # lower limit of 0.00001 (-5)

  # check if tidied
  if (!"factor" %in% names(F_matrix)) {
    cli::cli_abort(c(
      "{.var factor} column not detected:",
      "i" = "The F-matrix should contain a {.var factor} column.",
      "x" = "Did you forget to enable {.var tidy_output when reading the F-matrix}?"
    ))
  }

  # check if only one model_run
  if ("model_run" %in% names(F_matrix)) {
    if(length(unique(F_matrix$model_run)) > 1) {
      cli::cli_abort(c(
        "More than 1 run detected:",
        "i" = "The F-matrix contains more than 1 {.var model_run}.",
        "x" = "Did you select a base case run using a filter on {.var model_run}?"
      ))
    }
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
  if(!"factor" %in% class(F_matrix$species)) {
    F_matrix <- F_matrix %>%
      mutate(species = factor(paste0("`",species,"`")))
  }

  # create the plot.data
  df <- F_matrix %>%
    dplyr::mutate(value = dplyr::if_else(factor_profile == "concentration_of_species",
      dplyr::if_else(log10(value) < y_min, y_min, log10(value)),
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

  # prepare DISP data if needed

  # check if we need to do something with the DISP profiles so we prepare the
  # data here as the DISP avg needs to replace the concentration of species
  type_run_select <- "base_run"
  if(disp.errorbar) {

    if(!("DISP_avg" %in% unique(df$run_type))) {
      cli::cli_abort(c(
        "DISP results not found:",
        "i" = "The F-matrix does not contain a {.var run_type} column containing 'DISP_avg'.",
        "x" = "Did you provide the DISP results with the {.var F_matrix}?"
      ))
    }

    if(!("DISP_min" %in% unique(df$run_type))) {
      cli::cli_abort(c(
        "DISP results not found:",
        "i" = "The F-matrix does not contain a {.var run_type} column containing 'DISP_min'.",
        "x" = "Did you provide the DISP results with the {.var F_matrix}?"
      ))

    }

    if(!("DISP_max" %in% unique(df$run_type))) {
      cli::cli_abort(c(
        "DISP results not found:",
        "i" = "The F-matrix does not contain a {.var run_type} column containing 'DISP_max'.",
        "x" = "Did you provide the DISP results with the {.var F_matrix}?"
      ))


    }
    # set the correct profile for concentrations
    type_run_select <- "DISP_avg"
  }

  # order the x_labels
  if (!identical(xlabel.order, NA)) {
    df$species <- factor(as.character(df$species), levels = xlabel.order)
  }

  # create numeric column based on species. We use this to plot the data.
  df$numeric_x = as.numeric(df$species)

  # Calculate the upper axis value
  max_conc_log <- df %>%
    filter(factor_profile == "concentration_of_species",
           stringr::str_detect(run_type, type_run_select)) %>%
    select(value) %>%
    max()

  ## calculate the optimum axis
  breaks_log <- seq(y_min, ceiling(max_conc_log), 2)

  # if the breaks are to small, then increase them to the minimum
  if ((length(breaks_log) == 4) && (max(breaks_log = 1))) {
    breaks_log <- c(breaks_log, max(breaks_log) + 2)
  }

  labels_log <- unlist(lapply(breaks_log, function(x) {
    x.c <- paste0(ifelse(x >= 0, "1e+0", "1e-0"), ifelse(x >= 0, x, x * -1))
  }), use.names = FALSE)

  breaks_percentage <- seq(y_min, max(breaks_log), 1) * (length(seq(y_min, max(breaks_log), 1)) - 1) + y_min

  scale_items <- ceiling(length(breaks_percentage) / 2)

  sec.breaks <- c(seq(min(breaks_percentage), max(breaks_percentage), round(length(seq(min(breaks_percentage), max(breaks_percentage), 1)) / scale_items, 0)), max(breaks_percentage))
  sec.labels <- round(seq(0, 100, (100 / (length(sec.breaks) - 1))), 1)

  # scale the percentages to the secondary axis
  df <- df %>%
    dplyr::mutate(value = if_else(factor_profile == "percentage_of_species_sum", (value / 100) * (length(seq(y_min, max(breaks_log), 1)) - 1) + y_min, as.double(value)))

  plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = df %>%
                         filter(factor_profile == "concentration_of_species",
                                stringr::str_detect(run_type, type_run_select)),
                       ggplot2::aes(xmin = numeric_x - (bar.width/2),
                                    xmax = numeric_x + (bar.width/2),
                                    ymin = y_min, ymax = value),
                       col = bar.color,
                       fill = bar.color) +
    ggplot2::geom_point(ggplot2::aes(x = numeric_x, y = value),
                        data = df %>%
                          filter(factor_profile == "percentage_of_species_sum",
                                 stringr::str_detect(run_type, type_run_select)),
                        size = point.size,
                        col = point.color
    )

  # plot lollipop lines
  if (lollipops) {
    plot <- plot +
      ggplot2::geom_segment(
        data = df %>%
          filter(factor_profile == "percentage_of_species_sum",
                 stringr::str_detect(run_type, type_run_select)),
        ggplot2::aes(x = numeric_x,
                     xend = numeric_x,
                     y = y_min,
                     yend = value),
        color = point.color,
        alpha = 0.5
      )
  }

  # plot disp errorbars
  if(disp.errorbar) {
      disp.df <- df %>%
        filter(factor_profile == "concentration_of_species",
               stringr::str_detect(run_type, "DISP_")) %>%
        select(numeric_x, factor, run_type, value) %>%
      pivot_wider(names_from = run_type, values_from = value)

    plot <- plot +
      geom_errorbar(data = disp.df, aes(x = numeric_x,
                                        ymin=DISP_min,
                                        ymax=DISP_max),
                    width=0.4,
                    color = disp.errorbar.color,
                    size = 1,
                    position=position_dodge(.9))

  }

  plot <- plot +
    ggplot2::facet_grid(factor ~ .) +
    ggplot2::ylab(ylab) +
    ggplot2::scale_y_continuous(
      limits = c(y_min, max(breaks_log)),
      breaks = breaks_log,
      labels = labels_log,
      sec.axis = ggplot2::sec_axis(
        trans = ~ . * (length(seq(y_min, max(breaks_log), 1)) - 1) + y_min,
        name = ylab2,
        breaks = sec.breaks,
        labels = sec.labels
      ),
      expand = expansion(mult = c(0,0.05))
    ) +
    ggplot2::annotation_logticks(sides = "l") +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(name = xlab,
                              labels = parse(text = levels(df$species)),
                              breaks = sort(unique(df$numeric_x)),
                              expand = expansion(mult = expand.mult)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = xlabel.angle,
                                          size = x.font.size,
                                          hjust = xlabel.hjust,
                                          vjust = xlabel.vjust),
      axis.ticks.y.right = ggplot2::element_line(color = point.color),
      axis.line.y.right = ggplot2::element_line(color = point.color),
      axis.text.y.right = ggplot2::element_text(color = point.color),
      axis.title.y.right = ggplot2::element_text(color = point.color),
      panel.grid.major.x = element_blank()
    )

  # remove vertical grid lines?
  if (rm.grid.x) {
    plot <- plot +
      ggplot2::theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  }

  plot

  # should return more (class me2tools)
  #################
  # output
  #################
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
