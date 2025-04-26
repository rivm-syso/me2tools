#' Plot BS Boxplot using the EPA-PMF style
#'
#' This function plots a box and whisker plot based on the information read for
#' all BS runs using \code{me2_BS_read_F}. Note, that unlike the aggregated
#' results calculated by \code{me2_BS_read_F}, no filtering is applied by
#' default by this function. This means that each factor in each BS result is
#' used in the creation of the boxplot, even if their correlation is below
#' the defined correlation threshold. Factors with a correlation below the
#' threshold need to be removed manually based on the output provided by
#' \code{me2_BS_read_F}
#'
#' @param BS_results Tibble containing the F-matrix results for the BS runs.
#'   A list, containing data, as provided by \code{me2_BS_read_F} is also accepted.
#' @param xlabel.angle What angle should the x-axis labels be presented in? If
#'   your labels are long, \code{45} degrees can be useful, which is also the
#'   default.
#' @param xlabel.order The labels containing the species on the x-axis are
#'   plotted based on the factor levels. This parameter can contain the levels
#'   and the species are transformed to factors using these levels. N.B. these
#'   levels can contain code for expressions.
#' @param y_min The minimum value for the y-axis, which is used to cut-off the
#'   logarithmic scale. The default setting is \code{-5}, which corresponds to
#'   1xE-5.
#' @param xlab x-axis label.
#' @param ylab Primary y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
#' @param x.font.size The size of the xtick labels. Defaults to 10.
#' @param box.color Provide the fill color of the IQR box. Is set to
#'   \sQuote{burlywood1}, as an approximation of the color used in EPA-PMF
#'   program.
#' @param outlier.color The fill color for the outliers, with the
#'   default being \sQuote{firebrick}, similar to the color used in the EPA-PMF
#'   program.
#' @param outlier.size The size of the point of the outliers,
#'   defaults to 2.
#' @param outlier.shape The shape of the point of the outliers. The
#'   default is 3 (+).
#' @param base.color The fill color for the base results, with the
#'   default being \sQuote{royalblue}, similar to the color used in the EPA-PMF
#'   program.
#' @param base.size The size of the point for the base results,
#'   defaults to 2.
#' @param base.shape The shape of the point for the base results. The
#'   default is 23 (diamond).
#' @param show.plot A logical argument evaluating to TRUE or FALSE indicating
#'   whether the plot should be shown as default.
#' @param rm.grid.x Should the vertical grid lines be removed? In some cases
#'   there are a lot of species, causing the vertical grid lines to clutter the
#'   plot. If set to \code{TRUE} these lines are removed. Defaults to
#'   \code{FALSE}.
#' @param facet.parse.label Should the labels be parsed using the
#'   \code{labeller = label_parsed}? If set to \code{TRUE} then \code{"SO[2]"}
#'   will use subscript on the labels shown in the facet.
#' @param ... Other parameters, for example renamed parameters.
#'
#' @return me2tools list containing the ggplot2 with box plot results for the
#'   selected BS runs and the call to produce the plots.
#'
#' @export
#'
#' 
epa_plot_BSboxplot <- function(BS_results,
                               xlabel.angle = 45,
                               xlabel.order = NA,
                               y_min = -5,
                               ylab = "Concentration of species",
                               xlab = "Species",
                               auto.text = TRUE,
                               x.font.size = 10,
                               box.color = "burlywood1",
                               outlier.color = "firebrick",
                               outlier.size = 2,
                               outlier.shape = 3,
                               base.color = "royalblue",
                               base.size = 2,
                               base.shape = 23,
                               show.plot = TRUE,
                               rm.grid.x = FALSE,
                               facet.parse.label = FALSE,
                               ...) {
  # checks
  if (is.list(BS_results)) {
    if (exists("data", where = BS_results)) {
      BS.results <- BS_results$data
    } else {
      cli::cli_abort(c(
        "{.var data} not detected:",
        "i" = "The list {.var BS_results} should contain a {.var data} item.",
        "x" = "Did you provide the output of {.function me2_BS_read_F}?"
      ))
    }
  } else {
    BS.results <- BS_results
  }
  
  if(!"value" %in% names(BS.results)) {
    cli::cli_abort(c(
      "{.var value} not detected:",
      "i" = "The {.var data} in {.var BS_results} should contain a {.var value} column.",
      "x" = "Did you set {.var tidy_output} to {.var TRUE} when calling {.function me2_BS_read_F}?"
    ))
  }

  BS.base.results <- BS.results %>%
    filter(
      run_type == "BS_base",
      factor_profile == "concentration_of_species"
    )

  # prepare data
  BS.results <- BS.results %>%
    filter(
      run_type != "BS_base",
      factor_profile == "concentration_of_species"
    )

  if (nrow(BS.results) == 0) {
    cli::cli_abort(c(
      "{.var error estimates} not detected:",
      "i" = "The input does not contain {.var run_type} associated with the error estimates.",
      "x" = "Did you provide the output of {.function me2_BS_read_F}?"
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

  # order the x_labels
  if (!identical(xlabel.order, NA)) {
    BS.results$species <- factor(as.character(BS.results$species), 
                                 levels = xlabel.order)
  }
  if (!"factor" %in% class(BS.results$species)) {
    BS.results <- BS.results %>%
      mutate(species = factor(paste0("`", species, "`")))
  }

  x_labels <- levels(BS.results$species)

  plot.output <- ggplot2::ggplot(
    BS.results %>%
      dplyr::mutate(value = ifelse(value <= 10^y_min, 10^y_min, value)),
    aes(species, value),
    fill = "lightgray"
  ) +
    ggplot2::geom_boxplot(
      fill = box.color,
      outlier.color = outlier.color,
      outlier.shape = outlier.shape,
      outlier.size = outlier.size
    ) +
    ggplot2::geom_point(
      data = BS.base.results %>%
        mutate(value = ifelse(value <= 10^y_min, 10^y_min, value)),
      aes(x = species, y = value),
      color = base.color,
      fill = base.color,
      size = base.size,
      shape = base.shape
    )

  if (facet.parse.label) {
    plot.output <- plot.output +
      ggplot2::facet_grid(factor ~ ., labeller = label_parsed)
  } else {
    plot.output <- plot.output +
      ggplot2::facet_grid(factor ~ .)
  }

  plot.output <- plot.output +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(10^y_min, NA)
    ) +
    ggplot2::theme_bw() +
    annotation_logticks(sides = "l") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = xlabel.angle,
        hjust = 1,
        size = x.font.size
      ),
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(xlab,
      labels = parse(text = x_labels)
    ) +
    ylab(ylab)

  # remove vertical grid lines?
  if (rm.grid.x) {
    plot.output <- plot.output +
      ggplot2::theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  }

  # run garbage collector
  gc()


  ##################################################################
  ##                        Prepare output                        ##
  ##################################################################

  # print(metcor.plot)
  output <- list(
    "plot" = plot.output,
    call = match.call()
  )
  class(output) <- "me2tools"
  if (show.plot) {
    plot(output$plot)
  }
  invisible(output)
}
