#' MetCor plot options
#'
#' This function provides a default list of options supported by the
#' \code{metcor_plot} function.
#'
#' The \code{metcor_plot.options} contains five different list items for several
#' different parts of the plots. These items are \dQuote{plot}: main plot
#' settings; \dQuote{receptor}: settings for the receptor; \dQuote{raster}:
#' settings for the MetCor raster; \dQuote{legend}: settings for the legend;
#' \dQuote{annotation}: settings for the annotation.
#'
#' # Plot options
#'
#' List containing various options that are applicable for the total plot.
#'
#'  \itemize{
#'    \item \code{panel.background} is the color of the background of the plot.
#'      With countries this color is mostly applied to the sea. As such the
#'      default value for this color is \dQuote{#97dbf2}, but can be any HEX
#'      formatted code.
#'    \item \code{show.scale} A logical evaluating to TRUE or FALSE indicating
#'      whether a scale should be shown at the lower right hand side of the
#'      plot. Defaults to \code{TRUE}.
#'    \item \code{center.from} How to calculate the center point. Can be one of
#'      \dQuote{raster}, \dQuote{receptor} or \dQuote{manual}. In case of
#'      \dQuote{manual} the \code{center.point} has to be provided as well.
#'    \item \code{center.point} A vector containing lon and lat coordinates:
#'      \code{c("lon" = 0, "lat" = 0)} to be used as the center point.
#'    \item \code{zoom.level} This is the zoom level, works the same way as
#'      with OpenStreet map (OSM) in the sense that a larger number equals a
#'      larger zoom. Can be overridden by \code{xlim} and \code{ylim}.
#'    \item \code{xlim} limits of x, used for manual zoom. Defaults are set to
#'      \code{NA}.
#'    \item \code{ylim} limits of y, used for manual zoom. Defaults are set to
#'      \code{NA}.
#' }
#'
#' # Receptor options
#'
#' List containing various options that are applicable to plotting the receptor
#' spatial data frame (if provided).
#'
#' \itemize{
#'    \item \code{size} size of the receptor points in the map. Defaults to
#'      \code{3}.
#'    \item \code{shape} shape of the receptor points in the map. Defaults to
#'      \code{21},
#'        see the options for \code{ggplot2} for more information.
#'    \item \code{color} color of the stroke of the receptor points in the map.
#'        Defaults to \code{"blueviolet"}.
#'    \item \code{fill} color of the fill of the receptor points in the map.
#'        Defaults to \code{"blueviolet"}.
#'    \item \code{stroke} size of the stroke of the receptor points in the map.
#'        Defaults to \code{1}.
#' }
#'
#' # Raster options
#' List containing various options that are applicable to plotting the MetCor
#' raster results
#'
#' \itemize{
#'   \item \code{type} raster type, can be \code{discrete} or \code{gradient}.
#'     Defaults to \code{discrete}.
#'   \item \code{alpha} alpha transparency for the raster. Defaults to
#'     \code{0.65}.
#'   \item \code{smooth.factor} smoothing factor to linear smoothing. Defaults
#'     to \code{5}. When the value \code{0} is used, no smoothing is applied
#'   \item \code{gradient.color.palette} vector containing colors for the
#'     palette. Defaults to \code{pals::kovesi.rainbow(100)}, but basically
#'     every color palette can be used (e.g.
#'     \code{openair::openColours(scheme="jet")}).
#'   \item \code{discrete.breaks} legend breaks. Default setting is
#'     \code{automatic}, allowing the automatic calculation of equal bins based
#'     on the maximum of the raster. Can also be a vector containing hard coded
#'     breaks.
#'   \item \code{discrete.colors} vector containing color codes. Defaults to
#'     \code{automatic} which uses the following colors: \code{c("#2892c7",
#'     "#a0c79b", "#fafa64", "#fa8d34", "#e81014")}
#'   \item \code{layers.before} can contain a list with several simple feature
#'     (SF) layers that are plotted below the MetCor results. Each item in the
#'     list should contain a list with \dQuote{data}: dataframe with SF;
#'     \dQuote{color} color of the outline of the shapes; \dQuote{size} size of
#'     the outline of the shape; \dQuote{fill} the fill color of the shape.
#'   \item \code{layers.afer} can contain a list with several SF layers that
#'     are plotted on top of the MetCor results. See \code{layers.before} for
#'     requirements.
#' }
#'
#' # Legend options
#'
#' List containing various options that are applicable to the legend of the plot
#'
#' \itemize{
#'    \item \code{show} Should the legend be shown or hidden? Defaults to
#'      \code{TRUE}.
#'    \item \code{background.color} The color of the legend background.
#'    \item \code{background.alpha} The transparancy of the legend background
#'      color.
#'    \item \code{size} size of the legend. Defaults to \code{1.2}.
#'    \item \code{margins} List containing the margins for \dQuote{t,b,l,r}
#'      sides of the legend
#'    \item \code{title} title of the legend. Defaults to \code{NA}.
#'    \item \code{title.fontsize} font size of the legend. Defaults to \code{14}.
#'    \item \code{labels} vector containing the legend labels. If \code{NA}
#'      legend labels are automatically constructed using default setting
#'    \item \code{label.fontsize} font size of the labels. Defaults to
#'      \code{12}.
#'    \item \code{label.fontface} font face of the labels. Defaults to
#'      \code{plain}.
#'    \item \code{label.digits} number of digits shown in the labels. Defaults
#'      to \code{1}. Please note that this setting is only for the presentation
#'      of the labels, unrounded values are used in calculations
#'  }
#'
#' # Annotation options
#'
#' List containing various options that are applicable to the default plot
#' annotation option.
#'
#'  \itemize{
#'    \item \code{text} annotation text for the upper right corner.
#'    \item \code{fontsize} font size of the annotation. Defaults to \code{8}.
#'    \item \code{fontface} font face of the annotation. Defaults to \code{bold}.
#'  }
#'
#' @seealso \code{\link{metcor_import}}, \code{\link{metcor_project_raster}}, 
#' \code{\link{metcor_plot}}, \code{\link{metcor_export}}, 
#' \code{\link{metcor_fix_hysplit}}
#'
#' @export
#'
#' @importFrom pals kovesi.rainbow
#'


metcor_plot_options <- function() {
  ## Metcor.plot.options
  metcor.plot.options <- list(
    "plot" = list(
      "panel.background" = "#97dbf2",
      "show.scale" = TRUE,
      "center.from" = "raster",
      "center.point" = c("lon" = 0, "lat" = 0),
      "zoom.level" = 1,
      "xlim" = NA,
      "ylim" = NA
    ),
    "receptor" = list(
      "size" = 3,
      "shape" = 21,
      "color" = "blueviolet",
      "fill" = "blueviolet",
      "stroke" = 1
    ),
    "raster" = list(
      "type" = "discrete",
      "alpha" = 0.65,
      "smooth.factor" = 5, # use 0 for no smoothing (TRUE/FALSE remove!!)
      "gradient.color.palette" = pals::kovesi.rainbow(100), # renamed! grad.col.pal
      "discrete.breaks" = "automatic", # use min max or provide the breaks! renamed disc.breaks
      "discrete.colors" = "automatic", # use min max or provide the breaks! renamed disc.colors
      "layers.before" = list(),
      "layers.after" = list()
    ),
    "legend" = list(
      "show" = TRUE,
      "background.color" = "grey", # NEW!!
      "background.alpha" = 0.3, # NEW!!
      "size" = 1.2,
      "margins" = list(
        "t" = 7,
        "b" = 7,
        "r" = 7,
        "l" = 7,
        "unit" = "pt"
      ), # NEW!!
      "title" = NA,
      "title.fontsize" = 14,
      "title.fontface" = "bold", # NEW!!
      "labels" = NA,
      "label.fontsize" = 12,
      "label.fontface" = "plain",
      "label.digits" = 1
    ),
    "annotation" = list(
      "text" = NA,
      "fontsize" = 8,
      "fontface" = "bold"
    )
  )

  return(metcor.plot.options)
}
