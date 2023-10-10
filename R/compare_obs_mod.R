#' Compare the observed total mass vs the modelled contributions in the G-matrix
#'
#' Based on the observed total mass and the modelled contributions in the
#' G-matrix this function performs a simpe linear regression on the observations
#' and the sum of the factors. The function can perform an ordinary least 
#' squares regression (OLS) and a robust regression. In the latter regression 
#' outliers have less influence on the determination of the regression line. 
#' Also the regression can be based on a model with and without intercept. For 
#' usage of multilinear regression (MLR) see function 
#' \code{\link{compare_obs_mod_mlr}}
#'
#' @param data Tibble containing at least two columns of data: 1) the 
#'   measured mass concentrations (i.e., PM2.5 or PM10) and 2) the sum of the
#'   factors in concentration units. Additional information can be present, but
#'   will not be used.
#' @param x The name of the column in the "data" tibble containing the
#'   observations for the regression.
#' @param y The name of the column in the "data" tibble containing the
#'   sum of all factors for the regression.
#' @param robust Boolean used to switch between normal (OLS) or robust (RLS)
#'   regression.
#' @param maxit The limit on the number of IWLS iterations used when robust
#'   regression is performed. Defaults to 200. Can be increased in case of error
#'   messages related to failed convergence in the number of steps used.
#' @param intercept Should an intercept be used when doing the regression?
#'   Default is FALSE, forcing the regression through the origin.
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
#' @param geomtextpath This variable contains specific settings for the optional
#'   package geomtextpath. When plotting the mod.lines, the geomtextpath 
#'   package will add labels to the lines for clarification, but only if this
#'   package is installed. If it is not installed, the function will only plot
#'   the lines. The variable contains the labels and the hjust for positioning 
#'   the labels. The hjust shifts the label up and down the line, and is 
#'   entered as a value between 0-1. The labels are character labels that are
#'   plotted as is.
#'
#' @returns me2tools list containing the ggplot2 for the regression of the
#'  observations and the sum of the factors in concentration units (LR). The
#'  results of the LR models are also given, as well as the data used to 
#'  calculate these models. Additional information is provided regarding the 
#'  use of robust regression.
#'
#' @seealso \code{\link{compare_obs_mod_mlr}}
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
compare_obs_mod <- function(data,
                            x = "pm10",
                            y = "sum_factors",
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
                            asp.ratio.1 = FALSE,
                            geomtextpath = tibble("labels" = list("0.5" = "0.5",
                                                                  "1.0" = "1.0",
                                                                  "2.0" = "2.0"),
                                                  "hjust" = list("0.5" = 0.95,
                                                                 "1.0" = 0.95,
                                                                 "2.0" = 0.95))) {
  
  if (ncol(data) < 2) {
    cli::cli_abort(c(
      "Not enough information in {.var data}",
      "i" = "The {.var data} should have at least two columns.",
      "x" = "Did you forget the total mass or sum of factors?"
    ))
  }

  # check if x variable exists
  if (!x %in% names(data)) {
    cli::cli_abort(c(
      "{.var x} column not detected:",
      "i" = "The {.var data} should contain a {x} column.",
      "x" = "Did you provide the correct {.var x} to be used?"
    ))
  }
  
  if (!y %in% names(data)) {
    cli::cli_abort(c(
      "{.var y} column not detected:",
      "i" = "The {.var data} should contain a {y} column.",
      "x" = "Did you provide the correct {.var y} to be used?"
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

  # create the data for regression
  regr_data <- data %>% 
    rename(x := !!sym(x),
           y := !!sym(y))
  
  ##############################################################################
  # see if we can create a model using variables
  if (intercept) {
    fm <- y ~ x
  } else {
    fm <- y ~ x + 0
  }

  # calculate model
  if (robust) {
    modelLM <- MASS::rlm(fm, data = regr_data, maxit = maxit)
  } else {
    modelLM <- stats::lm(fm, data = regr_data)
  }
 
  # create plot
  n_data <- seq(min(regr_data$x), 
                max(regr_data$x), 
                length = ceiling(max(regr_data$x)-min(regr_data$x)) * 4)
  
  conf_int <- predict(modelLM, 
                      newdata = tibble::tibble(x = n_data),
                      interval='confidence', 
                      level=CI.level)
  pred_int <- predict(modelLM, 
                      newdata = tibble::tibble(x = n_data),
                      interval='prediction', 
                      level=PI.level)
  
  regression.out <- tibble(x = n_data,
                           fit = conf_int[,1],
                           conf_lwr = conf_int[,2],
                           conf_upr = conf_int[,3],
                           pred_lwr = pred_int[,2],
                           pred_upr = pred_int[,3])
  
  scatter <- ggplot2::ggplot(regr_data,
                             ggplot2::aes(x = x, y = y))

  # prediction intervals (PI)
  if (PI) {
    if (PI.type == "area") {
      scatter <- scatter + 
        geom_ribbon(data = regression.out,
                    aes(x = x,
                        y = pred_lwr,
                        ymin = pred_lwr, 
                        ymax = pred_upr), 
                    fill = PI.color,
                    alpha = PI.alpha)
      
    } else if (PI.type == "line") {
      scatter <- scatter + 
        geom_line(data = regression.out,
                  aes(x=x, y= pred_lwr), 
                  color = PI.color,
                  linetype = PI.linetype,
                  linewidth = PI.linewidth) +
        geom_line(data = regression.out,
                  aes(x=x, y= pred_upr), 
                  color = PI.color,
                  linetype = PI.linetype,
                  linewidth = PI.linewidth)
    } else {
      # wrong choice.
      cli::cli_abort(c(
        "{.var PI.type} not supported:",
        "i" = "The provided {.var PI.type} of {PI.type} is not supported.",
        "x" = "{.var PI.type} should be one of 'area' or 'line'."
      ))
    }
  }

  # add points and regression line
  scatter <- scatter +
    geom_point(data = regr_data,
               shape = point.shape,
               colour = point.color,
               alpha = point.alpha,
               size = point.size)
  
  
  # confidence intervals (CI)
  if (CI) {
    if (CI.type == "area") {
      scatter <- scatter + 
        geom_ribbon(data = regression.out,
                    aes(x = x,
                        y = conf_lwr,
                        ymin = conf_lwr, 
                        ymax = conf_upr), 
                    fill = CI.color,
                    alpha = CI.alpha)
      
    } else if (CI.type == "line") {
      
      scatter <- scatter + 
        geom_line(data = regression.out,
                  aes(x=x, y= conf_lwr), 
                  color = CI.color,
                  linetype = CI.linetype,
                  linewidth = CI.linewidth) +
        geom_line(data = regression.out,
                  aes(x=x, y= conf_upr), 
                  color = CI.color,
                  linetype = CI.linetype,
                  linewidth = CI.linewidth)
    } else {
      # wrong choice.
      cli::cli_abort(c(
        "{.var CI.type} not supported:",
        "i" = "The provided {.var CI.type} of {CI.type} is not supported.",
        "x" = "{.var CI.type} should be one of 'area' or 'line'."
      ))
    }
  }
  
  #add regression line
  scatter <- scatter +
      geom_line(data = regression.out,
                aes(x=x, y= fit), 
                color = line.color,
                linetype = "solid",
                linewidth = line.size)

  max_axis <- max(c(
    ceiling(ggplot_build(scatter)$layout$panel_scales_x[[1]]$range$range / 10) * 10,
    ceiling(ggplot_build(scatter)$layout$panel_scales_y[[1]]$range$range / 10) * 10
  ))
  
  min_axis <- min(c(
    ceiling(ggplot_build(scatter)$layout$panel_scales_x[[1]]$range$range / 10) * 10,
    ceiling(ggplot_build(scatter)$layout$panel_scales_y[[1]]$range$range / 10) * 10
  ))

  if (mod.line) {
    # check if package geomtextpath is installed
    if (system.file(package='geomtextpath')=="") {
      # not installed, so use default
      scatter <- scatter +
        ggplot2::geom_abline(intercept = 0,
                             slope = 1,
                             color = "black",
                             linetype = "solid",
                             linewidth = 0.25) +
        ggplot2::geom_abline(intercept = 0,
                             slope = 0.5,
                             color = "black",
                             linetype = "dashed",
                             linewidth = 0.25) +
        ggplot2::geom_abline(intercept = 0,
                             slope = 2,
                             color = "black",
                             linetype = "dashed",
                             linewidth = 0.25)
    } else {
      scatter <- scatter +
        geomtextpath::geom_labelabline(label = geomtextpath$labels[["0.5"]], 
                                       slope = 0.5, 
                                       intercept = 0, 
                                       hjust = geomtextpath$hjust[["0.5"]],
                                       boxcolour = NA,
                                       boxlinewidth = 0,
                                       fill = NA,
                                       gap = TRUE,
                                       linetype = "dashed",
                                       linewidth = 0.25) +
        geomtextpath::geom_labelabline(label = geomtextpath$labels[["1.0"]], 
                                       slope = 1, 
                                       intercept = 0, 
                                       hjust = geomtextpath$hjust[["1.0"]],
                                       boxcolour = NA,
                                       boxlinewidth = 0,
                                       fill = NA,
                                       gap = TRUE,
                                       linetype = "solid",
                                       linewidth = 0.25) +
        geomtextpath::geom_labelabline(label = geomtextpath$labels[["2.0"]], 
                                       slope = 2, 
                                       intercept = 0, 
                                       hjust = geomtextpath$hjust[["2.0"]],
                                       boxcolour = NA,
                                       boxlinewidth = 0,
                                       fill = NA,
                                       gap = TRUE,
                                       linetype = "dashed",
                                       linewidth = 0.25)
    }
  }

  if(regress.details) {
    if (robust) {
      scatter <- scatter +
        ggpmisc::stat_poly_eq(mapping = ggpmisc::use_label("eq"), 
                              formula = fm,
                              method = MASS::rlm,
                              method.args = list("maxit" = maxit))
    } else {
      scatter <- scatter +
        ggpmisc::stat_poly_eq(ggplot2::aes(label = paste("atop(", ggplot2::after_stat(eq.label), ",", ggplot2::after_stat(adj.rr.label), ")", sep = "")),
                              formula = fm,
                              method = "lm")
    }
  }
  
  scatter <- scatter +
    ggplot2::theme_bw() +
    xlab(xlab) +
    ylab(ylab)

  if (asp.ratio.1) {
    scatter <- scatter +
      ggplot2::coord_fixed(ratio = 1,
                           xlim = c(min_axis, max_axis),
                           ylim = c(min_axis, max_axis),
                           expand = FALSE,
                           clip = "on")
      
  }

  

  #################
  # output
  #################
  # print(metcor.plot)
  output <- list(
    "plot" = scatter,
    "formula_LR" = fm,
    "model_LR" = modelLM,
    "data_LR" = regr_data,
    "robust_regression" = robust,
    "regression.output" = regression.out,
    call = match.call()
  )
  class(output) <- "me2tools"
  if (show.plot) {
    plot(output$plot)
  }
  invisible(output)
}
