#' Plot contributions using the EPA-PMF style
#'
#' Function to plot the G_matrix contributions for a single run into faceted 
#' plots.
#'
#' @param G_matrix Tibble containing the G-matrix results for a single base run
#'   in normalised or concentration units.
#' @param unit Which units should be plotted? By default the units of the 
#'   G_matrix are \dQuote{normalised} to 1. The other option, which is also the
#'   default, is \dQuote{concentration}
#' @param xlabel.angle What angle should the x-axis labels be presented in? If
#'   your labels are long, \code{45} degrees can be useful, which is also the
#'   default.
#' @param xlabel.order The labels containing the species on the x-axis are
#'   plotted based on the factor levels. This parameter can contain the levels
#'   and the species are transformed to factors using these levels. N.B. these
#'   levels can contain code for expressions.
#' @param xlab x-axis label.
#' @param ylab Primary y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
#' @param x.font.size The size of the xtick labels. Defaults to 10.
#' @param show.plot A logical argument evaluating to TRUE or FALSE indicating
#'   whether the plot should be shown as default.
#' @param rm.grid.x Should the vertical grid lines be removed? In some cases
#'   there are a lot of species, causing the vertical grid lines to clutter the
#'   plot. If set to \code{TRUE} these lines are removed. Defaults to
#'   \code{FALSE}.
#' @param facet.parse.label Should the labels be parsed using the
#'   \code{labeller = label_parsed}? If set to \code{TRUE} then \code{"SO[2]"}
#'   will use subscript on the labels shown in the facet.
#'   
#' @return me2tools list containing the ggplot2 with contributions 
#'   plot and the call used to create this plot.
#'   

epa_plot_contributions <- function(G_matrix,
                                   unit = "concentration",
                                   xlabel.angle = 45,
                                   xlabel.order = NA,
                                   ylab = "Concentration of species",
                                   xlab = "Date",
                                   auto.text = TRUE,
                                   x.font.size = 10,
                                   show.plot = TRUE,
                                   rm.grid.x = FALSE,
                                   facet.parse.label = FALSE,
                                   ...) {
  
  ##################################################################
  ##                            Checks                            ##
  ##################################################################
  
  # check if tidied
  if (!"factor" %in% names(G_matrix)) {
    cli::cli_abort(c(
      "{.var factor} column not detected:",
      "i" = "The G-matrix should contain a {.var factor} column.",
      "x" = "Did you forget to enable {.var tidy_output} when reading the G-matrix?"
    ))
  }
  
  # check if only one model_run
  if ("model_run" %in% names(G_matrix)) {
    if (length(unique(G_matrix$model_run)) > 1) {
      cli::cli_abort(c(
        "More than 1 run detected:",
        "i" = "The G-matrix contains more than 1 {.var model_run}.",
        "x" = "Did you select a base case run using a filter on {.var model_run}?"
      ))
    }
  }
  
  G_matrix <- G_matrix  %>%
    filter(unit == unit)
  
  if (nrow(G_matrix) == 0) {
    cli::cli_abort(c(
      "{.var unit} does not have the right input:",
      "i" = "{.var unit} does not contain {unit}.",
      "x" = "Did you provide the correct input for the {.var G_matrix}?"
    ))
  }
  
  
  plot.output <- ggplot2::ggplot(
    data = G_matrix,
    ggplot2::aes(x = date, y = value)
  ) +
    ggplot2::geom_point(colour = "red", size = 1, alpha = 1)
  
  if (facet.parse.label) {
    plot.output <- plot.output +
      ggplot2::facet_grid(factor ~ ., scales = "free_y", labeller = label_parsed)  
    
  } else {
    plot.output <- plot.output +
      ggplot2::facet_grid(factor ~ ., scales = "free_y")  
  }
  
  plot.output <- plot.output +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xlabel.angle, 
                                              hjust = 1,
                                              size = x.font.size),
          legend.position="none",
          panel.grid.minor = element_blank()) +
    ylab(ylab) +
    xlab(xlab)
  
  
  
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