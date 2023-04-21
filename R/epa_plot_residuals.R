#' Function to plot histograms of the residuals with limit values.
#' 
#' This function plots the results of the scaled residuals using histograms.
#' 
#' @param me2_residuals The output of the function \code{me2_read_residuals}, 
#'   which reads the data file with residual results
#' @param limit.species A character vector containing the names of species that
#'   should be plotted. This is ideal to create bigger plots for examining the
#'   species. Defaults to \code{NA} so no limit is applied.
#' @param residual.limits The limits of the scaled residuals. These limits are
#'   plotted using vertical lines in the histogram.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
#' @param bin.number The number of bins used, defaults to 30,
#' @param bin.color The color of the bins, defaults to "steelblue",
#' @param bin.alpha The transparency of the bins, defaults to 0.7,
#' @param guide.color = The color of the guidelines based on the 
#'   residual.limits.Defaults to "red".
#' @param guide.linetype = The line type of the guidelines based on the 
#'   residual.limits.Defaults to "solid".
#' @param guide.linewidth The line width of the guidelines based on the 
#'   residual.limits.Defaults to 1.
#' @param x.font.size The size of the xtick labels. Defaults to 10.
#' @param xlabel.angle What angle should the x-axis labels be presented in? If
#'   your labels are long, \code{45} degrees can be useful, which is also the
#'   default.
#' @param facet.col The number of columns for the faceted plot. If the number of
#'   columns is set to 1, the faceting will be in rows only. Defaults to 3.
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
#' @return me2tools list containing the ggplot2 with scaled residual results 
#'   for the selected species and the call to produce the plot.   

epa_plot_residuals <- function(me2_residuals,
                               limit.species = NA,
                               residual.limits = c(-4,4),
                               xlab = "Scaled residuals",
                               ylab = "Number of residuals in group",
                               auto.text = TRUE,
                               bin.number = 30,
                               bin.color = "steelblue",
                               bin.alpha = 0.7,
                               guide.color = "red",
                               guide.linetype = "solid",
                               guide.linewidth = 1,
                               x.font.size = 10,
                               xlabel.angle = 45,
                               facet.col = 3,
                               show.plot = TRUE,
                               rm.grid.x = FALSE,
                               facet.parse.label = FALSE) { 

  
  
 
  # check if residual_scaled is present in the data
  if (!"residual_scaled" %in% unique(me2_residuals$residual_type)) {
    cli::cli_abort(c(
      "Scaled residuals not found:",
      "i" = "The {.var me2_residuals} should contain a {.var residual_type} with value 'residual_scaled'.",
      "x" = "Did you provide the correct {.var me2_residuals} to be used?"
    ))
  }
  
  # check if tidied data
  if (!"species" %in% names(me2_residuals)) {
    # tidy the data 
    
    species.order <- names(me2_residuals[4:length(me2_residuals)])

        me2_residuals <- me2_residuals %>% 
      pivot_longer(cols = -c("residual_type", "model_run", "date"),
                   names_to = "species",
                   values_to = "value") %>% 
      mutate(species = factor(species, levels = species.order))
    
  }
  
  # check the length of the model_run
  if (length(unique(me2_residuals$model_run)) > 1) {
    cli::cli_abort(c(
      "Multiple model runs detected:",
      "i" = "The {.var me2_residuals} should contain only one {.var model_run}.",
      "x" = "Did you forget to select a single 'base_run' from the {.var me2_residuals}?"
    ))
  }
  
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
  
  # filter on residual_scaled
  plot.data <- me2_residuals %>% 
    filter(residual_type == "residual_scaled")
  
  # generate list outside the limits
  outside.limits <- plot.data %>% 
    filter(!between(value, min(residual.limits), max(residual.limits)))
  
  # apply limit to species.
  if (length(limit.species) == 1) {
    if(!is.na(limit.species)) {
      plot.data <- plot.data %>% 
        filter(species %in% limit.species)
      
      if (nrow(plot.data) == 0) {
        cli::cli_abort(c(
          "No data to plot",
          "i" = "There is no data after applying {.var limit.species}.",
          "x" = "Did you provide the correct {.var limit.species} to be used?"
        ))
      }
    }
  } else {
    plot.data <- plot.data %>% 
      filter(species %in% limit.species)
    
    if (nrow(plot.data) == 0) {
      cli::cli_abort(c(
        "No data to plot",
        "i" = "There is no data after applying {.var limit.species}.",
        "x" = "Did you provide the correct {.var limit.species} to be used?"
      ))
    }
  }

  # create plot
  plot.output <- ggplot2::ggplot(data = plot.data,
                  aes(x=value)) +
  geom_histogram(alpha = bin.alpha, 
                 bins = bin.number, 
                 fill = bin.color) +
  geom_vline(xintercept=residual.limits, 
             colour = guide.color, 
             linetype = guide.linetype, 
             linewidth = guide.linewidth)
  
  if (facet.parse.label) {
    plot.output <- plot.output +
      ggplot2::facet_wrap(~species, ncol = facet.col, scales = "free_y", labeller = label_parsed)  
    
  } else {
    plot.output <- plot.output +
      ggplot2::facet_wrap(~species, ncol = facet.col, scales = "free_y")  
  }
  
  plot.output <- plot.output +
    theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xlabel.angle, 
                                                       hjust = 1,
                                                       size = x.font.size),
                   legend.position="none",
                   panel.grid.minor = element_blank())

  # generate output
  output <- list(
    "residual.plot" = plot.output,
    "outside.limits" = outside.limits,
    call = match.call()
  )
  class(output) <- "me2tools"
  
  if (show.plot) {
    plot(output$residual.plot)
  }
  
  invisible(output)
}