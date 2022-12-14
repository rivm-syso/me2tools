#' Create temporal plots of contributions (G-matrix)
#'
#' This function can create temporal plots for the factors in the G-contribution
#' matrix. Default uses all statistics as calculated by geom_boxplot.
#'
#' @param mydata A data frame of hourly (or higher temporal resolution data).
#'   Must include a \code{date} field and at least one variable (factor
#'   contributions) to plot.
#' @param factor Name of variable to plot, typically a factor contribution 
#'   from the G-matrix.
#' @param type \code{type} determines how the data are split i.e. conditioned,
#'   and then plotted. The \dQuote{default} will produce a single plot using the
#'   entire data. Type can be one of the built-in types as detailed in
#'   \code{cutData} e.g. \dQuote{season}, \dQuote{year}, \dQuote{weekday} and
#'   so on. For example, \code{type = "season"} will produce four plots --- one
#'   for each season.
#'
#'   It is also possible to choose \code{type} as another variable in the data
#'   frame. If that variable is numeric, then the data will be split into four
#'   quantiles (if possible) and labelled accordingly. If type is an existing
#'   character or factor variable, then those categories/levels will be used
#'   directly. This offers great flexibility for understanding the variation of
#'   different variables and how they depend on one another.
#' @param group This sets the grouping variable to be used. For example, if a
#'   data frame had a column \code{year} setting \code{group = "year"} will
#'   plot the years as individual box plots.
#' @param facet This sets the faceting variable to be used. For example, if a
#'   data frame had a column \code{site} setting \code{facet = "site"} will
#'   plot all sites individually in their own panels.
#' @param na.rm A logical evaluating to TRUE or FALSE indicating whether NA
#'   values should be stripped before the computation proceeds.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param auto.text Either \code{TRUE} (default) or \code{FALSE}. If
#'   \code{TRUE} titles and axis labels will automatically try and format
#'   pollutant names and units properly e.g.  by subscripting the \sQuote{2} in
#'   NO2.
#' @param txt.x.rm A logical evaluating to TRUE or FALSE indicating whether the
#'   x-axis tick labels should be removed.
#' @param txt.y.rm A logical evaluating to TRUE or FALSE indicating whether the
#'   y-axis tick labels should be removed.
#' @param cols Colours to be used for plotting. Options include
#'   \dQuote{default}, \dQuote{increment}, \dQuote{heat}, \dQuote{jet} and
#'   \code{RColorBrewer} colours --- see the \code{openair} \code{openColours}
#'   function for more details. For user defined the user can supply a list of
#'   colour names recognised by R (type \code{colours()} to see the full list).
#'   An example would be \code{cols = c("yellow", "green", "blue")}
#' @param whisk.lim User defined limits for the whiskers. Default is
#'   \code{NULL}, in this case the whiskers calculated by \dQuote{geom_boxplot}
#'   are used. In other cases a vector should be provided with two values
#'   between 0..1 denoting the percentile (i.e., \code{c(0.05, 0.95)} will use
#'   the P05 and P95 for the whisker limits).
#' @param point.size The size of the points. Defaults to 1.25.
#' @param line.size The size of the lines. Defaults to 1.
#' @param facet.col The number of columns for the faceted plot. If the number of
#'   columns is set to 1, the faceting will be in rows only. Defaults to 2.
#' @param facet.scales Should scales be fixed ("fixed", the default), free 
#'   ("free"), or free in one dimension ("free_x", "free_y")? See also 
#'   ggplot2::facet_wrap or ggplot2::facet_grid.
#' @param ... Other parameters passed onto \code{cutData}. For example, in the
#'   case of \code{cutData} the option \code{hemisphere = "southern"}.
#'
#' @return \code{temporal_factors} returns an object of class ``me2tools''.
#'   The object includes three main components: \code{call}, the command used
#'   to generate the plot; \code{plot}, the default plot based on the settings;
#'   and \code{box.plot}, the associated box and whisker plots in the style of
#'   Tukey.  If retained, e.g. using \code{output <- temporal_factors(mydata,
#'   "factor1")}, this output can be used to recover the data, reproduce or
#'   rework the original plot or undertake further analysis.
#'
#'   An me2tools output can be manipulated using a number of generic operations,
#'   including \code{print}, \code{plot} and \code{summary}.
#'
#' @seealso \code{\link[ggplot2]{geom_boxplot}}, \code{\link[openair]{cutData}}
#'   
#' @importFrom grDevices adjustcolor
#' @import cli
#' @import openair
#' @import grDevices
#' @import ggplot2
#' @import dplyr
#' @importFrom stats reformulate
#' @importFrom stats quantile
#' 
temporal_contributions <- function(mydata,
                                   factor = NA,
                                   type = "default",
                                   group = NA,
                                   facet = NA,
                                   na.rm = TRUE,
                                   xlab = NA,
                                   ylab = NA,
                                   auto.text = TRUE,
                                   txt.x.rm = FALSE,
                                   txt.y.rm = FALSE,
                                   cols = "hue",
                                   whisk.lim = NULL,
                                   point.size = 1.25,
                                   line.size = 1,
                                   facet.col = 2,
                                   facet.scales = "fixed",
                                   ...) {

  ## extra.args setup
  extra.args <- list(...)

  # check if type is length 1
  if (length(type) > 1) {
    cli::cli_abort(c(
      "Must have only one type:",
      "i" = "There {?is/are} {length(type)} type{?s}.",
      "x" = "You've tried to add multiple elements for {.var type}. Do you want
        to apply a facet or group variable instead?"
    ))
  }

  if((!identical(facet, NA)) & (!identical(group, NA))) {
    cli::cli_abort(c(
      "Cannot have both facets and grouping:",
      "i" = "Currently this function cannot handle facet and grouping variables.",
      "x" = "Both {.var facet} and {.var group} are provided."
    ))
  }

  # check if type exists in data if not part of cutData
  type_cutData = c("default",
                   "year",
                   "hour",
                   "month",
                   "season",
                   "weekday",
                   "weekend",
                   "monthyear",
                   "daylight",
                   "dst")

  if (!(type %in% type_cutData)) {
    # check if in dataframe
    if (!(type %in% names(mydata))) {
      cli_abort(c(
        "Must have a valid {.var type}:",
        "i" = "{.var type} does not exists in 'cutData' or in the data.",
        "x" = "Did you forget to add the column `{type}` to the data?"
      ))
    }
  }

  if(length(whisk.lim) > 0) {
    if(length(whisk.lim) == 1) {
      cli::cli_abort(c(
        "Must have two values:",
        "i" = "There {?is/are} {length(whisk.lim)} whisker limit{?s}.",
        "x" = "Two values for {.var whisk.lim} should be provided."
      ))
    }
    if(length(whisk.lim) > 2) {
      cli::cli_abort(c(
        "Must have only two values:",
        "i" = "There {?is/are} {length(whisk.lim)} whisker limit{?s}.",
        "x" = "You've tried to provide more that two values for {.var whisk.lim}."
      ))
    }
    if(max(whisk.lim) > 1) {
      cli::cli_abort(c(
        "Maximum cannot be larger than 1:",
        "i" = "The maximum for {.var whisk.lim} is {max(whisk.lim)}.",
        "x" = "Values for {.var whisk.lim} should be between [0..1]."
      ))
    }
    if(max(whisk.lim) < 0) {
      cli::cli_abort(c(
        "Minimum cannot be smaller than 0:",
        "i" = "The minimum for {.var whisk.lim} is {min(whisk.lim)}.",
        "x" = "Values for {.var whisk.lim} should be between [0..1]."
      ))
    }
  }

  # check if we have a faceting variable
  if (identical(facet, NA)) {
    if (identical(group, NA)) {
      vars <- type
    } else {
      vars <- c(group, type)
    }
    nfacet <- 1
  } else {
    if (!(facet %in% names(mydata))) {
      cli::cli_abort(c(
        "Must have a valid {.var facet}:",
        "i" = "{.var facet} does not exists in the data.",
        "x" = "Did you forget to add the column `{facet}` to the data?"
      ))
    }
    if (identical(group, NA)) {
      vars <- c(facet, type)
    } else {
      vars <- c(group, facet, type)
    }
    nfacet <- length(unique(mydata[[facet]])) ## number of facets
  }

  # create the colors
  myColors <- openair::openColours(cols, nfacet)
  boxplotColor <- grDevices::adjustcolor(myColors, alpha.f = 0.5)

  # cut the data
  mydata <- openair::cutData(mydata, type = vars)

  # Here is where we calculate the data for all the plots using the geom_boxplot
  # routine.
  # different plot routines for none faceting and faceting
  if (identical(facet, NA)) {
    # use geom_boxplot to create the default plot and grab the data.
    # we first create a white boxplot
    if (identical(group, NA)) {
      box.plot <- ggplot2::ggplot(
        mydata,
        ggplot2::aes(x = !!sym(type), y = !!sym(factor))
      ) +
        ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
      
      # then we grab the data
      box.data <- ggplot2::layer_data(box.plot)
  
      # since we have type variable, we need to apply this as well.
      # The output data has now a GROUP variable that we can use.
      box.data <- box.data %>%
        mutate(!!sym(type) := factor(group))
      # in order to make this work, we need to set the correct group names
      # to the panel.
      if ("factor" %in% class(mydata[[type]])) {
        levels(box.data[[type]]) <- levels(mydata[[type]])
      } else {
        # not a factor
        levels(box.data[[type]]) <- sort(unique(mydata[[type]]))
      }
  
      ############################################
      # see if we need to replace the whiskers
      ############################################
      if(length(whisk.lim) == 2) {
        additional.data <- mydata %>%
          dplyr::select(!!vars, !!factor) %>%
          dplyr::group_by(dplyr::across(all_of(vars))) %>%
          dplyr::summarise(
            n_wllimit = stats::quantile(!!sym(factor), probs = min(whisk.lim), na.rm = na.rm),
            n_wulimit = stats::quantile(!!sym(factor), probs = max(whisk.lim), na.rm = na.rm),
            .groups = "drop_last"
          )
        # combine with the boxplot data
        box.data <- box.data %>%
          left_join(., additional.data, by = vars) %>%
          rename(ymin_old = ymin,
                 ymax_old = ymax,
                 ymin = n_wllimit,
                 ymax = n_wulimit)
  
        cli::cli_warn(c(
          "The limits of the whiskers have been changed by the user.",
          "i" = "Whisker limits are set to {whisk.lim}, which are probabilities for quantiles."
        ))
      }
  
      # now we can construct the plots.
      default.plot <- ggplot2::ggplot(
        data = box.data,
        ggplot2::aes(x = !!sym(type), y = middle, group = 1)
      ) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), fill = myColors, colour = NA, alpha = 0.25) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = myColors, colour = NA, alpha = 0.25) +
        ggplot2::geom_line(colour = myColors, size = line.size) +
        ggplot2::geom_point(colour = myColors, size = point.size) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
  
      # make sure
      #test.data <- mydata %>%
      #  dplyr::select(!!vars, !!factor)
      #default.plot + geom_boxplot(mapping = ggplot2::aes(x = !!sym(type), y = !!sym(factor), group = !!sym(type), alpha = 0.2), outlier.shape = NA, data = test.data)
  
      # boxplot (replacing intial boxplot used to get data)
      box.plot <- ggplot2::ggplot(
        box.data,
        ggplot2::aes(x = !!sym(type),
            ymin = ymin,
            lower = lower,
            middle = middle,
            upper = upper,
            ymax = ymax)) +
        geom_boxplot(stat = "identity", color = "white", outlier.shape = NA, na.rm = TRUE) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
  
      # create the colors
      iqrColor <- adjustcolor(myColors, alpha.f = 0.5)
      # add colored segments to plot
      box.plot <- box.plot +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = lower, xend = xmax, yend = lower),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = lower, xend = xmin, yend = upper),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmax, y = lower, xend = xmax, yend = upper),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = upper, xend = xmax, yend = upper),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = middle, xend = xmax, yend = middle),
          color = myColors, size = 1.2
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = x, y = upper, xend = x, yend = ymax),
          color = "black"
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = x, y = ymin, xend = x, yend = lower),
          color = "black"
        )
    } else { # check grouping
      # we have to apply grouping here.
      box.plot <- ggplot2::ggplot(
        mydata,
        ggplot2::aes(x = !!sym(type), y = !!sym(factor), fill = !!sym(group))
      ) +
        ggplot2::geom_boxplot(outlier.shape = NA, na.rm = TRUE) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
      
      # then we grab the data
      box.data <- ggplot2::layer_data(box.plot)
      
      # since we have type variable, we need to apply this as well.
      # The output data has now a GROUP variable that we can use.
      box.data <- box.data %>%
        mutate(!!sym(type) := factor(group),
               !!sym(group) := factor(group))
      # in order to make this work, we need to set the correct group names
      # to the panel. Since we have grouped data we need to find the number
      # of items for each type.

      type_items <- mydata %>% 
        select(!!sym(type), !!sym(group)) %>% 
        unique() %>% 
        group_by(!!sym(type)) %>% 
        tally()
      
      types <- c() # emtpy to start with
      groups <- c() # emtpy to start with
      
      # also grab the group_levels
      if ("factor" %in% class(mydata[[group]])) {
        group_levels <- levels(mydata[[group]])
      } else {
        # not a factor
        group_levels <- sort(unique(mydata[[group]]))
      }

      for (type_level in levels(mydata[[type]])) {
        group_count <- as.numeric(type_items %>% 
                                    filter(!!sym(type) == type_level) %>% 
                                    select(n))
        types <- c(types, replicate(group_count, type_level))
        # add groups
        for (item in seq(1,group_count,1)) {
          groups <- c(groups, group_levels[item])
        }
      }

      if ("factor" %in% class(mydata[[type]])) {
        box.data[[type]] <- factor(types, levels(mydata[[type]]))
      } else {
        # not a factor
        box.data[[type]] <- factor(types, sort(unique(mydata[[type]])))
      }
      if ("factor" %in% class(mydata[[group]])) {
        box.data[[group]] <- factor(groups, levels(mydata[[group]]))
      } else {
        # not a factor
        box.data[[groups]] <- factor(groups, sort(unique(mydata[[groups]])))
      }

      ############################################
      # see if we need to replace the whiskers
      ############################################
      if(length(whisk.lim) == 2) {
        additional.data <- mydata %>%
          dplyr::select(!!vars, !!factor) %>%
          dplyr::group_by(dplyr::across(all_of(vars))) %>%
          dplyr::summarise(
            n_wllimit = stats::quantile(!!sym(factor), probs = min(whisk.lim), na.rm = na.rm),
            n_wulimit = stats::quantile(!!sym(factor), probs = max(whisk.lim), na.rm = na.rm),
            .groups = "drop_last"
          )
        # combine with the boxplot data
        box.data <- box.data %>%
          left_join(., additional.data, by = vars) %>%
          rename(ymin_old = ymin,
                 ymax_old = ymax,
                 ymin = n_wllimit,
                 ymax = n_wulimit)
        
        cli::cli_warn(c(
          "The limits of the whiskers have been changed by the user.",
          "i" = "Whisker limits are set to {whisk.lim}, which are probabilities for quantiles."
        ))
      }
      # Change MyColors to the default fill colors
      myColors <- box.data$fill
      
      # now we can construct the plots.
      #default.plot <- ggplot2::ggplot(
      #  data = box.data,
      #  ggplot2::aes(x = !!sym(type), y = middle)
      #) +
      #  ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax, fill = myColors, group = !!sym(group)), colour = NA, alpha = 0.25) +
      #  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = myColors, group = !!sym(group)), colour = NA, alpha = 0.25) +
      #  ggplot2::geom_line(ggplot2::aes(group = !!sym(group)),colour = myColors, size = line.size) +
      #  ggplot2::geom_point(ggplot2::aes(group = !!sym(group)),colour = myColors, size = point.size) +
      #  ggplot2::theme_bw() +
      #  ggplot2::theme(legend.position = "none")
      default.plot = NULL
      
      # Change MyColors to the default fill colors
      myColors <- box.data$fill
      
      # boxplot (replacing intial boxplot used to get data)
      box.plot <- ggplot2::ggplot(
        box.data,
        ggplot2::aes(x = !!sym(type),
                     ymin = ymin,
                     lower = lower,
                     middle = middle,
                     upper = upper,
                     ymax = ymax,
                     fill = !!sym(group))) +
        geom_boxplot(stat = "identity", color = NA, outlier.shape = NA, na.rm = TRUE) +
        scale_fill_manual(values = replicate(length(group_levels), NA), 
                          na.value = NA) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")

      # create the colors
      iqrColor <- adjustcolor(myColors, alpha.f = 0.5)
      # add colored segments to plot
      box.plot <- box.plot +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = lower, xend = xmax, yend = lower, group = !!sym(group)),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = lower, xend = xmin, yend = upper, group = !!sym(group)),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmax, y = lower, xend = xmax, yend = upper, group = !!sym(group)),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = upper, xend = xmax, yend = upper, group = !!sym(group)),
          color = iqrColor
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = xmin, y = middle, xend = xmax, yend = middle, group = !!sym(group)),
          color = myColors, size = 1.2
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = x, y = upper, xend = x, yend = ymax, group = !!sym(group)),
          color = "black"
        ) +
        ggplot2::geom_segment(
          data = box.data,
          mapping = ggplot2::aes(x = x, y = ymin, xend = x, yend = lower, group = !!sym(group)),
          color = "black"
        )
      
    } # end grouping

    # end
  } else { # faceted data
    # we first create a white boxplot
    box.plot <- ggplot2::ggplot(
      mydata,
      ggplot2::aes(x = !!sym(type), y = !!sym(factor))
    ) +
      ggplot2::geom_boxplot(color = "black", outlier.shape = NA, na.rm = TRUE) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::facet_wrap(stats::reformulate(facet), ncol = facet.col)

    # then we grab the data
    box.data <- ggplot2::layer_data(box.plot)
    # since we have group variable, we need to apply this as well.
    # The output data has now a PANEL variable that we can use.
    box.data <- box.data %>%
      mutate(!!sym(type) := factor(group),
             !!sym(facet) := factor(PANEL))
    # in order to make this work, we need to set the correct type and facet names
    # to the panel.
    if ("factor" %in% class(mydata[[type]])) {
      levels(box.data[[type]]) <- levels(mydata[[type]])
    } else {
      # not a factor
      levels(box.data[[type]]) <- sort(unique(mydata[[type]]))
    }
    if ("factor" %in% class(mydata[[facet]])) {
      levels(box.data[[facet]]) <- levels(mydata[[facet]])
    } else {
      # not a factor
      levels(box.data[[facet]]) <- sort(unique(mydata[[facet]]))
    }

    ############################################
    # see if we need to replace the whiskers
    ############################################
    if(length(whisk.lim) == 2) {
      additional.data <- mydata %>%
        dplyr::select(!!vars, !!factor) %>%
        dplyr::group_by(dplyr::across(all_of(vars))) %>%
        dplyr::summarise(
          n_wllimit = stats::quantile(!!sym(factor), probs = min(whisk.lim), na.rm = na.rm),
          n_wulimit = stats::quantile(!!sym(factor), probs = max(whisk.lim), na.rm = na.rm),
          .groups = "drop_last"
        )
      # combine with the boxplot data
      box.data <- box.data %>%
        left_join(., additional.data, by = vars) %>%
        rename(ymin_old = ymin,
               ymax_old = ymax,
               ymin = n_wllimit,
               ymax = n_wulimit)

      cli::cli_warn(c(
        "The limits of the whiskers have been changed by the user.",
        "i" = "Whisker limits are set to {whisk.lim}, which are probabilities for quantiles."
      ))
    }

    # now we can construct the plots.
    default.plot <- ggplot2::ggplot(
      data = box.data,
      ggplot2::aes(x = !!sym(type), y = middle, group = 1)
    ) +
      ggplot2::scale_fill_manual(values = myColors) +
      ggplot2::scale_color_manual(values = myColors) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = ymax), alpha = 0.25) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.25) +
      ggplot2::geom_line(ggplot2::aes(colour = !!sym(facet)), size = line.size) +
      ggplot2::geom_point(ggplot2::aes(colour = !!sym(facet)), size = point.size) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::aes(fill = !!sym(facet))

    # apply wrapping
    if (facet.col == 1) {
      default.plot <- default.plot +
        ggplot2::facet_grid(rows = vars(!!sym(facet)), scales = facet.scales)
    } else {
      default.plot <- default.plot +
        ggplot2::facet_wrap(stats::reformulate(facet), ncol = facet.col, , scales = facet.scales)
    }
    
    # add boxplot
    # we first create a white boxplot
    box.plot <- ggplot2::ggplot(
      box.data,
      ggplot2::aes(x = !!sym(type),
                   ymin = ymin,
                   lower = lower,
                   middle = middle,
                   upper = upper,
                   ymax = ymax,
                   group = group)) +
      geom_boxplot(stat = "identity", color = NA) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
    
    # apply wrapping
    if (facet.col == 1) {
      box.plot <- box.plot +
        ggplot2::facet_grid(rows = vars(!!sym(facet)), scales = facet.scales)
    } else {
      box.plot <- box.plot +
        ggplot2::facet_wrap(stats::reformulate(facet), ncol = facet.col, scales = facet.scales)
    }

    # create the colors
    iqrColor <- adjustcolor(myColors, alpha.f = 0.5)

    iqrColor <- rep(iqrColor[1], each = nrow(box.data))

    medianColor <- rep(myColors[1], each = nrow(box.data))

    # add colored segments to plot
    box.plot <- box.plot +
      ggplot2::geom_segment(
        data = box.data,
        mapping = ggplot2::aes(x = xmin, y = lower, xend = xmax, yend = lower),
        color = iqrColor
      ) +
      ggplot2::geom_segment(
        data = box.data,
        mapping = ggplot2::aes(x = xmin, y = lower, xend = xmin, yend = upper),
        color = iqrColor
      ) +
      ggplot2::geom_segment(
        data = box.data,
        mapping = ggplot2::aes(x = xmax, y = lower, xend = xmax, yend = upper),
        color = iqrColor
      ) +
      ggplot2::geom_segment(
        data = box.data,
        mapping = ggplot2::aes(x = xmin, y = upper, xend = xmax, yend = upper),
        color = iqrColor
      ) +
      ggplot2::geom_segment(
        data = box.data,
        mapping = ggplot2::aes(x = xmin, y = middle, xend = xmax, yend = middle),
        color = medianColor, size = 1.2
      ) +
      ggplot2::geom_segment(
        data = box.data,
        mapping = ggplot2::aes(x = x, y = upper, xend = x, yend = ymax),
        color = "black"
      ) +
      ggplot2::geom_segment(
        data = box.data,
        mapping = ggplot2::aes(x = x, y = ymin, xend = x, yend = lower),
        color = "black"
      )
  }
  if (!identical(xlab, NA)) {
    if (auto.text) {
      default.plot <- default.plot +
        ggplot2::xlab(openair::quickText(xlab))
      box.plot <- box.plot +
        ggplot2::xlab(openair::quickText(xlab))
    } else {
      if (!identical(xlab, NULL)) {
        default.plot <- default.plot +
          ggplot2::xlab(xlab)
        box.plot <- box.plot +
          ggplot2::xlab(xlab)
      } else {
        # remove the label
        default.plot <- default.plot +
          theme(axis.title.x = element_blank())
        box.plot <- box.plot +
          theme(axis.title.x = element_blank())
      }
    }
  }
  if (!identical(ylab, NA)) {
    if (auto.text) {
      default.plot <- default.plot +
        ggplot2::ylab(openair::quickText(ylab))
      box.plot <- box.plot +
        ggplot2::ylab(openair::quickText(ylab))
    } else {
      if (!identical(ylab, NULL)) {
        default.plot <- default.plot +
          ggplot2::ylab(ylab)
        box.plot <- box.plot +
          ggplot2::ylab(ylab)
      } else {
        # remove the label
        default.plot <- default.plot +
          theme(axis.title.y = element_blank())
        box.plot <- box.plot +
          theme(axis.title.y = element_blank())
      }
    }
  }
  # fix the hour levels when facet is set
  if ((!identical(facet, NA)) & (type == "hour")) {
    default.plot <- default.plot +
      ggplot2::scale_x_discrete(labels = c(
        "00", "", "02", "", "04", "",
        "06", "", "08", "", "10", "",
        "12", "", "14", "", "16", "",
        "18", "", "20", "", "22", ""
      ))

    box.plot <- box.plot +
      ggplot2::scale_x_discrete(labels = c(
        "00", "", "02", "", "04", "",
        "06", "", "08", "", "10", "",
        "12", "", "14", "", "16", "",
        "18", "", "20", "", "22", ""
      ))
  }

  if (txt.x.rm) {
    default.plot <- default.plot +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
    box.plot <- box.plot +
      ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }

  if (txt.y.rm) {
    default.plot <- default.plot +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
    box.plot <- box.plot +
      ggplot2::theme(axis.text.y = ggplot2::element_blank())
  }

  #  default.plot <- default.plot +
  #  #theme(axis.text.x = element_text(angle = 45, hjust = 0.5))
  #    ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(n.dodge = 2))


  output <- list(
    "plot" = default.plot,
    "box.plot" = box.plot,
    "data" = mydata,
    "box.data" = box.data,
    call = match.call()
  )
  class(output) <- "me2tools"
  invisible(output)
}
