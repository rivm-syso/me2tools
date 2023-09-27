#' Overview of the input of ME2
#' 
#' This function uses the input for ME2 to provide an overview of the species 
#' concentrations and uncertainties.
#' 
#' @param me2_input The output of the function \code{me2_read_input}, which
#'   reads an input file defined with a date-time column, followed by 
#'   alternating concentrations and uncertainties for a number of unknown 
#'   species. At the very least this variable is a list, containing a 
#'   \code{conc} and an \code{unc} dataframe.
#' @param concentration The name of the list variable containing the 
#'   concentrations. Defaults to \dQuote{conc}.
#' @param uncertainty The name of the list variable containing the 
#'   uncertainties. Defaults to \dQuote{unc}.
#' @param limit.species A character vector containing the names of species that
#'   should be plotted. This is ideal to create bigger plots for examining the
#'   species. Defaults to \code{NA} so no limit is applied.
#' @param sn.colors This parameter is a list containing the colors used in the
#'   concentration vs uncertainty plots. In these plots the calculated S/N ratio
#'   is used to color the points
#' @param sn.limits These are the limits for the S/N ratio. Values above the 
#'   \code{u_limit} are classified as "strong", the values below the
#'   \code{l_limit} are classified as "bad" and the values between the limits
#'   are classified as "weak". The default limits are \code{u_limit = 1} and
#'   \code{l_limit = 0.5} which correspond to the S/N calculation introduced in
#'   EPA-PMF version 5. Note this is only a classification and does not mean 
#'   that data should be excluded or included by default.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
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
#' @return me2tools list containing the ggplot2 with concentration vs 
#'   uncertainty for the selected species and the call to produce the plot.   


input_overview  <- function(me2_input,
                            concentration = "conc",
                            uncertainty = "unc",
                            limit.species = NA,
                            sn.colors = list("bad" = "firebrick3",
                                             "weak" = "goldenrod2",
                                             "strong" = "green4"),
                            sn.limits = list("u_limit" = 1,
                                             "l_limit" = 0.5),
                            xlab = "Concentration (ug/m3)",
                            ylab = "Uncertainty (ug/m3)",
                            auto.text = TRUE,
                            x.font.size = 10,
                            xlabel.angle = 45,
                            facet.col = 3,
                            show.plot = TRUE,
                            rm.grid.x = FALSE,
                            facet.parse.label = FALSE){
  
  # write checks
  
  # check if concentration and uncertainty exists.
  if (!concentration %in% names(me2_input)) {
    cli::cli_abort(c(
      "{.var me2_input} does not contain concentrations.",
      "i" = "The {.var concentration} = '{concentration}' variable could not be found in {.var me2_input}.",
      "x" = "Did you provide the correct {.var concentration} to be used?"
    ))
  }
  
  if (!uncertainty %in% names(me2_input)) {
    cli::cli_abort(c(
      "{.var me2_input} does not contain uncertainties",
      "i" = "The {.var uncertainty} = '{uncertainty}' variable could not be found in {.var me2_input}.",
      "x" = "Did you provide the correct {.var uncertainty} to be used?"
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

  species.order <- names(me2_input[[concentration]][2:length(me2_input[[concentration]])])
  
  conc.data <- me2_input[[concentration]] %>% 
    tidyr::pivot_longer(cols = -"date",
                        names_to = "species",
                        values_to = "value") %>% 
    mutate(type = "pmf.concentration",
           species = factor(species, levels = species.order))
  
  unc.data <- me2_input[[uncertainty]] %>% 
    tidyr::pivot_longer(cols = -"date",
                        names_to = "species",
                        values_to = "value") %>% 
    mutate(type = "pmf.uncertainty",
           species = factor(species, levels = species.order))
  
  plot.data <- dplyr::bind_rows(conc.data,
                                unc.data) %>% 
    pivot_wider(id_cols = c("date", "species"),
                names_from = "type",
                values_from = "value")
  
  epa.summary <- plot.data %>% 
    group_by(species) %>% 
    summarize(sn = epa_sn(x=pmf.concentration, x_unc = pmf.uncertainty),
              min = min(pmf.concentration, na.rm = TRUE),
              P25 = epa_percentile(pmf.concentration, prob = 0.25),
              P50 = epa_percentile(pmf.concentration, prob = 0.50),
              P75 = epa_percentile(pmf.concentration, prob = 0.75),
              max = max(pmf.concentration, na.rm = TRUE),
    ) %>% 
    mutate(guidance = dplyr::if_else(sn > sn.limits$u_limit, 
                                     "strong", 
                                     dplyr::if_else(sn < sn.limits$l_limit, 
                                                    "bad", 
                                                    "weak")))
  
  plot.data <- left_join(plot.data,
                         epa.summary,
                         by = join_by(species)) %>% 
    mutate(guidance = factor(guidance, levels = c("bad", "weak", "strong")))
  
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

  plot.output <- ggplot(data = plot.data, aes(x=pmf.concentration, y=pmf.uncertainty)) +
    geom_point(aes(color = guidance), size = 1) +
    scale_color_manual(values=c("bad" = sn.colors[["bad"]], 
                                "weak" = sn.colors[["weak"]], 
                                "strong" = sn.colors[["strong"]]))
  
  if (facet.parse.label) {
    plot.output <- plot.output +
      ggplot2::facet_wrap(~species, ncol = facet.col, scales = "free", labeller = label_parsed)  
    
  } else {
    plot.output <- plot.output +
      ggplot2::facet_wrap(~species, ncol = facet.col, scales = "free")  
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
    "cu.plot" = plot.output,
    "stats.table" = epa.summary,
    call = match.call()
  )
  class(output) <- "me2tools"
  
  if (show.plot) {
    plot(output$cu.plot)
  }
  
  invisible(output)
}