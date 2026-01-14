#' Plot Q/Qexp plot using the EPA-PMF style
#'
#' Using this function the residuals of the PMF solution are plotted in a
#' Q/Qexp plot (both profiles and contributions). This is a great way to
#' examine the residuals, specifically the species or samples that were not well
#' modelled (i.e., Q/Qexp > 2).
#'
#' @param residuals Tibble from `me2_read_residuals()` containing
#'   the results for the scaled residuals. In order to use this function only
#'   the results from a single base run should be provided.
#' @param num_factors Compulsory number of factors used to generate the
#'   residuals.
#' @param num_strong_species Number of species that are designated as "strong"
#'   variables in the analysis. Defaults to \code{NA}, meaning all species are
#'   set to "strong".
#' @param threshold The threshold that should be used to check the Q/Qexp
#'   against. Defaults to 2, as mentioned in the EPA-PMF manual as being the
#'   threshold to observe species or samples that have not been well modelled.
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
#' @param ylab Primary y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
#' @param date.breaks When this parameter is a string it is passed to the
#'   function \code{scales::breaks_width()} and deals with the number of x
#'   ticks. It contains the distance between each break. Either a number, or for
#'   date/times, a single string of the form "\{n\} \{unit\}", e.g. "1 month",
#'   "5 days". Unit can be of one "sec", "min", "hour", "day", "week",
#'   "month", "year". It can also be a single integer, in which case the
#'   variable is passed to \code{scales::breaks_pretty()}.
#' @param x.font.size The size of the xtick labels. Defaults to 10.
#' @param x.n.dodge The number of rows on the x-axis that should be used to
#'   render the labels. Useful for labels that otherwise overlap, for example
#'   with a large number of species. Also allows for a larger
#'   \sQuote{x.font.size}.
#' @param bar.color Provide the fill color for the Q/Qexp based bars for each
#'   species. Is set to \sQuote{cadetblue3}, as an approximation
#'   of the color used in EPA-PMF program.
#' @param bar.width The width of the bar, expressed as a value between \[0,1\].
#' @param bar.alpha The alpha transparency of the fill color of the bar,
#'   expressed as a value between \[0,1\].
#' @param point.color The color for the Q/Qexp point of the samples, with the
#'   default being \sQuote{red}.
#' @param point.size The size for the Q/Qexp point of the samples,
#'   defaults to 1.
#' @param point.alpha The alpha for the Q/Qexp point of the samples,
#'   defaults to 1.
#' @param threshold.color The color used for plotting the threshold as a
#'   horizontal line.
#' @param threshold.size The line size of the horizontal threshold line.
#'   Defaults to 1.
#' @param threshold.type The line type of the horizontal threshold line.
#'   Defaults to 2 (dashed line).
#' @param show.plot A logical argument evaluating to TRUE or FALSE indicating
#'   whether the plot should be shown by default.
#' @param show.legend A logical argument evaluating to TRUE or FALSE indicating
#'   whether the legend should be shown by default.
#' @param rm.grid.x Should the vertical grid lines be removed? In some cases
#'   there are a lot of species, causing the vertical grid lines to clutter the
#'   plot. If set to \code{TRUE} these lines are removed. Defaults to
#'   \code{FALSE}.
#' @param ... Other parameters, for example renamed parameters.
#'
#' @return me2tools list containing the ggplot2 with segment and point
#'   geometries and the call.
#'
#' @export
#'
#'
epa_plot_QQexp <- function(residuals,
                           num_factors = NA,
                           num_strong_species = NA,
                           threshold = 2,
                           xlabel.angle = 45,
                           xlabel.vjust = 1,
                           xlabel.hjust = 1,
                           xlabel.order = NA,
                           ylab = "Q/Qexp",
                           auto.text = TRUE,
                           date.breaks = NA,
                           x.font.size = 10,
                           x.n.dodge = 1,
                           bar.color = "cadetblue3",
                           bar.width = 0.9,
                           bar.alpha = 0.8,
                           point.color = "red",
                           point.size = 1,
                           point.alpha = 1,
                           threshold.color = "black",
                           threshold.size = 1,
                           threshold.type = 2,
                           show.plot = TRUE,
                           show.legend = FALSE,
                           rm.grid.x = FALSE,
                           ...) {
  # checks
  if (identical(num_factors, NA)) {
    cli::cli_abort(c(
      "{.var num_factors} not provided:",
      "i" = "The number of factors used to generate the residual data needs to be provided.",
      "x" = "{.var num_factors} should be larger than 1."
    ))
  }

  scaled_residuals <- residuals %>%
    filter(residual_type == "residual_scaled")

  if (nrow(scaled_residuals) == 0) {
    cli::cli_abort(c(
      "{.var residuals} needs to have 'residual_scaled",
      "i" = "{.var residuals} does not contain {.var residual_type} = 'residual_scaled'.",
      "x" = "{.var residual_type} should contain 'residual_scaled'."
    ))
  }

  # variables needed
  num_species <- ncol(scaled_residuals) - 3
  num_samples <- nrow(scaled_residuals)

  if (identical(num_strong_species, NA)) {
    num_strong_species <- num_species
  }


  ##################################################################
  ##                        Calculate Qexp                        ##
  ##################################################################
  Qexp <- (num_samples * num_strong_species) - ((num_factors * num_samples) + (num_factors * num_strong_species))

  ##################################################################
  ##                   Process scaled residuals                   ##
  ##################################################################

  # perhaps easier to have tidied data, as you only need to square one column
  scaled_residuals[, 4:ncol(scaled_residuals)] <- scaled_residuals[, 4:ncol(scaled_residuals)]^2

  ##################################################################
  ##                     Create F Q/Qexp plot                     ##
  ##################################################################
  Q.Qexp_F <- tibble(
    "species" = names(scaled_residuals[, 4:ncol(scaled_residuals)]),
    "Q.Qexp" = colSums(scaled_residuals[, 4:ncol(scaled_residuals)]) / (Qexp / num_strong_species)
  )

  # order the x_labels
  if (!identical(xlabel.order, NA)) {
    Q.Qexp_F$species <- factor(as.character(Q.Qexp_F$species), levels = xlabel.order)
  } else {
    Q.Qexp_F$species <- factor(Q.Qexp_F$species, levels = Q.Qexp_F$species)
  }

  plot_F <- ggplot2::ggplot(data = Q.Qexp_F, aes(x = species, y = Q.Qexp)) +
    ggplot2::geom_bar(
      stat = "identity",
      col = bar.color,
      fill = bar.color,
      width = bar.width,
      alpha = bar.alpha
    ) +
    ggplot2::geom_hline(
      yintercept = threshold,
      linetype = threshold.type,
      color = threshold.color,
      linewidth = threshold.size
    ) +
    ggplot2::scale_x_discrete(
      guide = ggplot2::guide_axis(angle = xlabel.angle, n.dodge = x.n.dodge)
    ) +
    ggplot2::scale_y_continuous(
      expand = expansion(mult = c(0, 0.05))
    ) +
    ylab(ylab) +
    xlab(NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "right",
      legend.margin = margin(b = -10),
      axis.text.x = ggplot2::element_text(
        size = x.font.size,
        hjust = xlabel.hjust,
        vjust = xlabel.vjust
      )
    )


  ##################################################################
  ##                     Create G Q/Qexp plot                     ##
  ##################################################################
  Q.Qexp_G <- tibble(
    "date" = scaled_residuals$date,
    "Q.Qexp" = rowSums(scaled_residuals[, 4:ncol(scaled_residuals)]) / (Qexp / num_samples)
  )

  plot_G <- ggplot2::ggplot(Q.Qexp_G, aes(x = date, y = Q.Qexp)) +
    ggplot2::geom_point(
      colour = point.color,
      size = point.size,
      alpha = point.alpha
    ) +
    ggplot2::geom_hline(
      yintercept = threshold,
      linetype = threshold.type,
      color = threshold.color,
      linewidth = threshold.size
    ) +
    ylab(ylab) +
    xlab(NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "right",
      legend.margin = margin(b = -10)
    )


  # remove legend?
  if (!show.legend) {
    plot_F <- plot_F +
      ggplot2::theme(legend.position = "none")
    plot_G <- plot_G +
      ggplot2::theme(legend.position = "none")
  }

  # remove vertical grid lines?
  if (rm.grid.x) {
    plot_F <- plot_F +
      ggplot2::theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
    plot_G <- plot_G +
      ggplot2::theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  }

  if (!identical(date.breaks, NA)) {
    if ("numeric" %in% class(date.breaks)) {
      plot_G <- plot_G +
        scale_x_datetime(breaks = scales::breaks_pretty(date.breaks))
    } else {
      plot_G <- plot_G +
        scale_x_datetime(breaks = scales::breaks_width(date.breaks))
    }
  }

  # find the species > threshold
  Q.Qexp_F_threshold <- Q.Qexp_F %>%
    filter(Q.Qexp > threshold)

  Q.Qexp_G_threshold <- Q.Qexp_G %>%
    filter(Q.Qexp > threshold)


  # run garbage collector
  gc()
  # if (!identical(dev.list(), NULL)) {
  #  dev.off()
  # }

  ##################################################################
  ##                        Prepare output                        ##
  ##################################################################

  # print(metcor.plot)
  output <- list(
    "Q.Qexp_F_plot" = plot_F,
    "Q.Qexp_G_plot" = plot_G,
    "Q.Qexp_F" = Q.Qexp_F,
    "Q.Qexp_G" = Q.Qexp_G,
    "Q.Qexp_F_threshold" = Q.Qexp_F_threshold,
    "Q.Qexp_G_threshold" = Q.Qexp_G_threshold,
    call = match.call()
  )
  class(output) <- "me2tools"
  if (show.plot) {
    plot(output$Q.Qexp_F_plot)
  }
  invisible(output)
}
