#' Plot error summary using the EPA-PMF style
#'
#' Function to plot the error summary using a combined F_matrix containing the
#' base run and at least one or all of the error estimates BS, DISP or BSDISP.
#' The results for the error estimates need to have a run_type of 
#' \dQuote{DISP_min} and \dQuote{DISP_max} for DISP; the \dQuote{BS_P05} and 
#' \dQuote{BS_P95} for BS and \dQuote{BSDISP_P05} and \dQuote{BSDISP_P95} for 
#' BSDISP, respectively. For the average base_run concentrations the 
#' concentration of each species from the base run is used.
#'
#' @param F_matrix Tibble containing the F-matrix results for the base runs 
#'   runs and at least one of the error estimates BS, DISP or BSDISP.
#' @param xlabel.angle What angle should the x-axis labels be presented in? If
#'   your labels are long, \code{45} degrees can be useful, which is also the
#'   default.
#' @param xlabel.parse should the xlabels be parsed as an expression? If set to
#'   \code{FALSE}, the labels are plotted as characters.
#' @param xlabel.vjust vertical justification of the xlabel, between \[0,1\].
#'   The default is NA, so that ggplot uses some heuristics to pick the best
#'   value for this parameter. Any other value is processed by changing theme
#'   settings. Note that for this to work the \code(xlabel.hjust) also needs to
#'   have a value.
#' @param xlabel.hjust horizontal justification of the xlabel, between \[0,1\].
#'   The default is \code{NA}, so that ggplot uses some heuristics to pick the
#'   best value for this parameter. Any other value is processed by changing
#'   theme settings. Note that for this to work the \code(xlabel.vjust) also 
#'   needs to have a value.
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
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} attempts will be made to automatically format titles,
#'   axis labels, pollutant names and units properly, e.g., by subscripting
#'   the \sQuote{2} in NO2.
#' @param x.font.size The size of the xtick labels. Defaults to 10.
#' @param x.n.dodge The number of rows on the x-axis that should be used to 
#'   render the labels. Useful for labels that otherwise overlap, for example
#'   with a large number of species. Also allows for a larger 
#'   \sQuote{x.font.size}.
#' @param bar.colors The bar colors for each error estimate (BS, DISP, BSDISP)
#'   and the color for the average of the base run.
#' @param bar.width The width of the bar, expressed as a value between \[0,1\].
#' @param show.plot A logical argument evaluating to TRUE or FALSE indicating
#'   whether the plot should be shown by default.
#' @param show.legend A logical argument evaluating to TRUE or FALSE indicating
#'   whether the legend should be shown by default.  
#' @param expand.mult Vector of multiplicative range expansion factors used on
#'   the x-axis . Defaults to \code{c(0.015,0.005)}, as this seems to work best.
#'   If length is 1, then both left and right x-axis are multiplied with this
#'   factor. If length is 2, the left limit is expanded by \code{expand.mult[1]}
#'   and the right by \code{expand.mult[2]}.
#' @param rm.grid.x Should the vertical grid lines be removed? In some cases
#'   there are a lot of species, causing the vertical grid lines to clutter the
#'   plot. If set to \code{TRUE} these lines are removed. Defaults to
#'   \code{FALSE}.
#' @param facet.parse.label Should the labels be parsed using the
#'   \code{labeller = label_parsed}? If set to \code{TRUE} then \code{"SO[2]"}
#'   will use subscript on the labels shown in the facet.
#' @param ... Other parameters, for example renamed parameters.
#'   
#' @return me2tools list containing the ggplot2 with error estimates summary 
#'   plot and the call used to create this plot.
#'   
#' @export
#'
#' @import cli
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' 

epa_plot_errorsummary <- function(F_matrix,
                             xlabel.angle = 45,
                             xlabel.parse = FALSE,
                             xlabel.vjust = NA,
                             xlabel.hjust = NA,
                             xlabel.order = NA,
                             xlabel.top.bottom = FALSE,
                             y_min = -5,
                             ylab = "Concentration of species",
                             xlab = "Species",
                             auto.text = TRUE,
                             x.font.size = 10,
                             x.n.dodge = 1,
                             bar.colors = list("BS" = "goldenrod2",
                                               "BSDISP" = "green4",
                                               "DISP" = "royalblue2",
                                               "base" = "black"),
                             bar.width = 0.9,
                             show.plot = TRUE,
                             show.legend = FALSE,
                             expand.mult = c(0.015, 0.005),
                             rm.grid.x = FALSE,
                             facet.parse.label = FALSE,
                             ...) {
  
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

  # check if it contains "concentration_of_species" and
  # "percentage_of_species_sum"
  
  # check if species is factor
  if (!"factor" %in% class(F_matrix$species)) {
    F_matrix <- F_matrix %>%
      mutate(species = factor(paste0("`", species, "`")))
  }
  
  
  ##################################################################
  ##                         Extract data                         ##
  ##################################################################
  
  id_variables <- c("factor_profile", "model_run", "species", "model_type", "group_type", "factor")
  
  df <- F_matrix %>% 
    dplyr::filter(factor_profile == "concentration_of_species") %>% 
    tidyr::separate(run_type, c("group_type", "plot_type"), "_") %>% 
    dplyr::mutate(plot_type = str_replace(plot_type, "run", "mean")) %>% 
    dplyr::mutate(plot_type = str_replace(plot_type, "avg", "mean")) %>% 
    dplyr::mutate(plot_type = str_replace(plot_type, "median", "mean")) %>% 
    dplyr::mutate(plot_type = str_replace(plot_type, "P95", "max")) %>% 
    dplyr::mutate(plot_type = str_replace(plot_type, "P05", "min")) %>% 
    tidyr::pivot_wider(id_cols = dplyr::all_of(id_variables),
                       names_from = "plot_type",
                       values_from = "value")
  
  
  #################################################################
  ##             Plot preparations (axis breaks etc)             ##
  #################################################################
  
  # order the x_labels
  if (!identical(xlabel.order, NA)) {
    df$species <- factor(as.character(df$species), levels = xlabel.order)
  }
  
  # create numeric column based on species. We use this to plot the data.
  df$numeric_x <- as.numeric(df$species)
  
  # get the segment data
  segments <- df %>%
    filter(group_type != "base")
  
  if (nrow(segments)==0) {
    stop("Did you provide error estimates?")
  }
  
  segments <- segments %>%
    mutate(mean = ifelse(mean <= 10^y_min, 10^y_min, mean),
           min = ifelse(min <= 10^y_min, 10^y_min, min),
           max = ifelse(max <= 10^y_min, 10^y_min, max),
           group_type = factor(group_type, levels = c("BS", "BSDISP", "DISP"))) %>% 
    droplevels(.) # drop unused levels
  
  
  # calculate the postions
  group.types <- levels(segments$group_type)
  num.segments <- length(group.types)
  
  if (num.segments == 0) {
    stop("nothing to plot!")
  }
  individual.bar.width <- bar.width/num.segments
  
  # define plot
  plot.output <- ggplot2::ggplot()
  plot.data <- tibble()
  
  group.type <- group.types[1]
  index <- 1
  for (group.type in group.types) {
    # with three segments, we need to correct the mid point
    tmp.plot.data <- segments %>% 
      filter(group_type == group.type)
    
    if (num.segments == 3) {
      if (index == 1) {
        tmp.plot.data <- tmp.plot.data %>% 
          dplyr::mutate(xmin = numeric_x - (individual.bar.width + (individual.bar.width/2)),
                        xmax = numeric_x - (individual.bar.width/2))
      } else if(index == 2) {
        tmp.plot.data <- tmp.plot.data %>% 
          dplyr::mutate(xmin = numeric_x - (individual.bar.width/2),
                        xmax = numeric_x + (individual.bar.width/2))
      } else if(index == 3) {
        tmp.plot.data <- tmp.plot.data %>% 
          dplyr::mutate(xmin = numeric_x + (individual.bar.width/2),
                        xmax = numeric_x + (individual.bar.width + (individual.bar.width/2)))
      }
    } else if (num.segments == 2) {
      if (index == 1) {
        tmp.plot.data <- tmp.plot.data %>% 
          dplyr::mutate(xmin = numeric_x - (individual.bar.width),
                        xmax = numeric_x)
      } else if(index == 2) {
        tmp.plot.data <- tmp.plot.data %>% 
          dplyr::mutate(xmin = numeric_x,
                        xmax = numeric_x + (individual.bar.width))
      }
    } else {
      if (index == 1) {
        tmp.plot.data <- tmp.plot.data %>% 
          dplyr::mutate(xmin = numeric_x - (individual.bar.width/2),
                        xmax = numeric_x + (individual.bar.width/2))
      }
    }
    
    if (nrow(plot.data) == 0) {
      plot.data <- tmp.plot.data 
    } else {
      plot.data <- dplyr::bind_rows(plot.data,
                                    tmp.plot.data)
    }
    index <- index + 1
  }

  plot.output <- plot.output +
    ggplot2::geom_rect(
      data = plot.data,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = min, 
        ymax = max,
        fill = group_type,
      ),
      #fill = bar.colors[[group.type]]
    )
  
  # add the base line
  base.line <- df %>%
    filter(group_type == "base") %>% 
    mutate(mean = ifelse(mean <= 10^y_min, 10^y_min, mean),
           group_type = factor(group_type, 
                               levels = c("BS", "BSDISP", "DISP"))) %>% 
    droplevels(.) # drop unused levels
  
  
  plot.output <- plot.output +
    geom_errorbar(
      data = base.line, aes(
        x = numeric_x,
        ymin = mean,
        ymax = mean,
        color = "base"
      ),
      width = bar.width,
      #color = bar.colors[["base"]],
      linewidth = 1,
      position = position_dodge(.9)
    )
  
  if (facet.parse.label) {
    plot.output <- plot.output +
      ggplot2::facet_grid(factor ~ ., labeller = label_parsed)  
    
  } else {
    plot.output <- plot.output +
      ggplot2::facet_grid(factor ~ .)  
  }
  
  # get the labels
  x_labels <- levels(df$species)
  
  if (xlabel.parse) {
    x_labels <- parse(text = x_labels)
  }
  
  # split into even and odd
  x_labels_top <- x_labels
  x_labels_bottom <- x_labels
  # define quick functions
  odd <- function(x)
    x %% 2 != 0
  even <- function(x)
    x %% 2 == 0
  # split the labels
  x_labels_top[odd(seq(1, length(x_labels_top)))] <- NA
  x_labels_bottom[even(seq(1, length(x_labels_bottom)))] <- NA
  
  
  plot.output <- plot.output +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(10^y_min,NA)
    ) +
    ggplot2::theme_bw() +
    annotation_logticks(sides = 'l') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = xlabel.angle, 
                                                       hjust = 1,
                                                       size = x.font.size),
          panel.grid.minor = element_blank())
  
  if ((identical(xlabel.vjust, NA)) &
      (identical(xlabel.hjust, NA))) {
    if (xlabel.top.bottom) {
      plot.output <- plot.output +
        ggplot2::scale_x_continuous(
          name = xlab,
          #labels = parse(text = x_labels_top),
          labels = x_labels_top,
          breaks = sort(unique(base.line$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = guide_axis(angle = xlabel.angle, 
                             n.dodge = x.n.dodge),
          sec.axis = ggplot2::sec_axis(
            trans = ~ .,
            breaks = sort(unique(base.line$numeric_x)),,
            #labels =  parse(text = x_labels_bottom),
            labels =  x_labels_bottom,
            guide = ggplot2::guide_axis(angle = xlabel.angle, 
                                        n.dodge = x.n.dodge)
          )
        ) +
        theme(axis.text.x.top = ggplot2::element_text(hjust = 0,
                                                      vjust = 0.5))
    } else {
      plot.output <- plot.output +
        ggplot2::scale_x_continuous(
          name = xlab,
          #labels = parse(text = x_labels),
          labels = x_labels,
          breaks = sort(unique(base.line$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = ggplot2::guide_axis(angle = xlabel.angle, 
                                      n.dodge = x.n.dodge)
        )
    }
  } else {
    if (xlabel.top.bottom) {
      plot.output <- plot.output +
        ggplot2::scale_x_continuous(
          name = xlab,
          #labels = parse(text = x_labels_top),
          labels = x_labels_top,
          breaks = sort(unique(base.line$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = guide_axis(angle = xlabel.angle, 
                             n.dodge = x.n.dodge),
          sec.axis = ggplot2::sec_axis(
            trans = ~ .,
            breaks = sort(unique(base.line$numeric_x)),,
            #labels =  parse(text = x_labels_bottom),
            labels =  x_labels_bottom,
            guide = ggplot2::guide_axis(n.dodge = x.n.dodge)
          )
        ) +
        theme(axis.text.x.top = ggplot2::element_text(hjust = xlabel.hjust,
                                                      vjust = xlabel.vjust))
    } else {
      plot.output <- plot.output +
        ggplot2::scale_x_continuous(
          name = xlab,
          #labels = parse(text = x_labels),
          labels = x_labels,
          breaks = sort(unique(base.line$numeric_x)),
          expand = expansion(mult = expand.mult),
          guide = ggplot2::guide_axis(n.dodge = x.n.dodge)
        )
    }
    plot.output <- plot.output +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = xlabel.angle,
          size = x.font.size,
          hjust = xlabel.hjust,
          vjust = xlabel.vjust
        ))
  }

  plot.output <- plot.output +
    ylab(ylab)  +
    scale_fill_manual(
      labels = c("BS" = "BS",
                 "BSDISP" = "BSDISP",
                 "DISP" = "DISP",
                 "base" = "BASE"),
      values = c("BS" = bar.colors[["BS"]], 
                 "BSDISP" = bar.colors[["BSDISP"]],
                 "DISP" = bar.colors[["DISP"]],
                 "base" = bar.colors[["base"]]),
      breaks = c("BS", 
                 "BSDISP",
                 "DISP",
                 "base")) +
    scale_color_manual(
      labels = c("BS" = "BS",
                 "BSDISP" = "BSDISP",
                 "DISP" = "DISP",
                 "base" = "BASE"),
      values = c("BS" = bar.colors[["BS"]], 
                 "BSDISP" = bar.colors[["BSDISP"]],
                 "DISP" = bar.colors[["DISP"]],
                 "base" = bar.colors[["base"]]),
      breaks = c("BS", 
                 "BSDISP",
                 "DISP",
                 "base")) +
      guides(color = guide_legend(order = 1,
                                  title = NULL),
             fill = guide_legend(order = 2,
                                  title = NULL))+
    theme(legend.position = "top",
          legend.justification = "right",
          legend.margin=margin(b=-10)) 
  
  # remove legend?
  if(!show.legend) {
    plot <- plot +
      ggplot2::theme(legend.position = "none")
  }

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