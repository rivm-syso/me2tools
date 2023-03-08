#' Plot the MetCor results
#'
#' This function plots the MetCor results. For the world map data from the
#' \code{rnaturalearth} package is used.
#'
#' @param metcor.raster MetCor raster object. Can be the default raster after
#'   reading (see \link{metcor_import}) or a transformed projected raster (see
#'   \link{metcor_project_raster})
#' @param metcor.plot.options List containing all the plot options (see
#'   \link{metcor_plot_options})
#' @param receptor Tibble containing the location of the receptors with
#'   longitude and latitude coordinates. Please note that the column names
#'   should be "lat" and "lon".
#' @param verbose show detailed messages of what is happening. Defaults to
#'   \code{TRUE}
#' @param show.plot A logical evaluating to TRUE or FALSE indicating whether the
#'   plot should be outputed as default.
#'
#' @return \code{metcor_plot} returns an object of class ``me2tools''.
#'   The object includes five main components: \code{call}, the command used
#'   to generate the plot; \code{plot}, the default plot based on the settings;
#'   \code{plot.options}, containing all information needed to generate the
#'   plot, \code{proj_crs}, the used CRS and \code{origin}, the center point
#'   used in the plot. If retained, e.g. using
#'   \code{output <- metcor_plot(metcor.raster)}, this output can be used to
#'   recover the data, reproduce or rework the original plot or undertake
#'   further analysis.
#'
#'   An me2tools output can be manipulated using a number of generic operations,
#'   including \code{print}, \code{plot} and \code{summary}.
#'
#'  @seealso \code{\link{metcor_import}}, \code{\link{metcor_project_raster}},
#'  \code{\link{metcor_plot_options}}
#'
#' @export
#'
#' @seealso \code{\link{metcor_import}}, \code{\link{metcor_project_raster}},
#' \code{\link{metcor_plot_options}}, \code{\link{metcor_export}},
#' \code{\link{metcor_fix_hysplit}}
#'
#' @importFrom classInt classIntervals
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 coord_sf
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_sf
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 unit
#' @importFrom ggspatial annotation_scale
#' @importFrom sf st_as_sf
#' @importFrom sf st_coordinates
#' @importFrom sf st_point
#' @importFrom sf st_read
#' @importFrom sf st_sfc
#' @importFrom sf st_transform
#' @importFrom sf st_make_grid
#' @importFrom sf st_set_crs
#' @importFrom stats na.omit
#' @importFrom stats density
#' @importFrom terra classify
#' @importFrom terra disagg
#' @importFrom terra ext
#' @importFrom terra res
#' @importFrom terra trim
#' @importFrom terra minmax
#' @importFrom tidyterra geom_spatraster
#' @importFrom pals kovesi.rainbow
#' @importFrom raster extent
#' @importFrom ggtext element_markdown
#' @importFrom ggnewscale new_scale_color
#' @import rnaturalearth
#' @import rnaturalearthdata
#' @import cli
#'
metcor_plot <- function(metcor.raster,
                        metcor.plot.options = metcor_plot_options(),
                        receptor = NA,
                        verbose = TRUE,
                        show.plot = TRUE) {

  if (!metcor.plot.options$raster$type %in% c("discrete", "gradient")){
    cli::cli_abort(c(
      "Invalid option for {.var raster$type} in the plot options:",
      "x" = "Only 'gradient' or 'discrete' is supported."
    ))
  }
  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                               OPTIONS                               ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  if (verbose) {
    message("- processing options")
  }
  ##################################################################
  ##                        Raster options                        ##
  ##################################################################
  if (metcor.plot.options$raster$type == "discrete") {
    raster.minmax <- terra::minmax(metcor.raster)
    if (length(metcor.plot.options$raster$discrete.breaks) == 1) {
      if (metcor.plot.options$raster$discrete.breaks == "automatic") {
        # calculate based on the max of the raster
        metcor.plot.options$raster$discrete.breaks <- c(
          raster.minmax[1, 1],
          raster.minmax[2, 1] / 16,
          raster.minmax[2, 1] / 8,
          raster.minmax[2, 1] / 4,
          raster.minmax[2, 1] / 2,
          raster.minmax[2, 1]
        )
      } else {
        cli::cli_abort(c(
          "Invalid option for {.var discrete.breaks}:",
          "x" = "Only 'automatic' or a vector with breaks are supported."
        ))
      }
    }
    if (length(metcor.plot.options$raster$discrete.colors) == 1) {
      if (metcor.plot.options$raster$discrete.colors == "automatic") {
        metcor.plot.options$raster$discrete.colors <- c(
          "#2892c7",
          "#a0c79b",
          "#fafa64",
          "#fa8d34",
          "#e81014"
        )
      } else {
        cli::cli_abort(c(
          "Invalid option for {.var discrete.colors}:",
          "x" = "Only 'automatic' or a vector with colors are supported."
        ))
      }
    }
    ## Make sure that the max value in the disc breaks matches the max in the
    ## raster
    if (max(metcor.plot.options$raster$discrete.breaks) < raster.minmax[2, 1]) {
      cli::cli_warn(c(
        "Maximum in data exceeds maximum `discrete.breaks` value:",
        "i" = "The highest value in the data ({raster.minmax[2, 1]}) is larger
          than the maximum value in the legend
          ({max(metcor.plot.options$raster$discrete.breaks)})."
      ))
    }
    if (min(metcor.plot.options$raster$discrete.breaks) > raster.minmax[1, 1]) {
      cli::cli_warn(c(
        "Minimum in data exceeds minimum `discrete.breaks` value:",
        "i" = "The lowest value in the data ({raster.minmax[1, 1]}) is smaller
          than the minimum value in the legend
          ({min(metcor.plot.options$raster$discrete.breaks)})."
      ))
    }
  }
  ## Calculate stats
  stats <- c(
    terra::global(metcor.raster, min, na.rm = TRUE)[1, 1],
    terra::global(metcor.raster, mean, na.rm = TRUE)[1, 1],
    terra::global(metcor.raster, max, na.rm = TRUE)[1, 1]
  )

  names(stats) <- c("min", "mean", "max")

  #################################################################
  ##                    Setting theme options                    ##
  #################################################################
  if (metcor.plot.options$raster$type == "discrete") {
    # different theme options for discrete legends
    theme_opts <- list(ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = metcor.plot.options$plot$panel.background, colour = NA
      ),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        size = 1
      ),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = ggplot2::unit(22, "points")
      ),
      legend.justification = c(0, 0),
      legend.position = c(0, 0),
      legend.title = ggtext::element_markdown(
        size = metcor.plot.options$legend$title.fontsize,
        lineheight = 1
      ),
      # legend.key = ggplot2::element_rect(colour= "black", size = 1),
      legend.text = ggtext::element_markdown(
        size = ggplot2::unit(metcor.plot.options$legend$label.fontsize, "points"),
        face = metcor.plot.options$legend$label.fontface
      ),
      legend.key.size = ggplot2::unit(
        metcor.plot.options$legend$size, "lines"
      ),
      legend.margin = ggplot2::margin(
        t = metcor.plot.options$legend$margins$t,
        l = metcor.plot.options$legend$margins$l,
        b = metcor.plot.options$legend$margins$b,
        r = metcor.plot.options$legend$margins$r,
        unit = metcor.plot.options$legend$margins$unit
      ),
      legend.background = ggplot2::element_rect(
        fill = ggplot2::alpha(
          metcor.plot.options$legend$background.color,
          metcor.plot.options$legend$background.alpha
        ),
        size = 0, linetype = 0
      )
    ))
  } else{
    # different theme options for gradient legends
    theme_opts <- list(ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = metcor.plot.options$plot$panel.background, colour = NA
      ),
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        size = 1
      ),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = ggplot2::unit(22, "points")
      ),
      legend.justification = c(0, 0),
      legend.position = c(0, 0),
      legend.title = ggtext::element_markdown(
        size = metcor.plot.options$legend$title.fontsize,
        lineheight = 1
      ),
      # legend.key = ggplot2::element_rect(colour= "black", size = 1),
      #legend.text = ggtext::element_markdown(
      legend.text = ggplot2::element_text(
        size = ggplot2::unit(metcor.plot.options$legend$label.fontsize, "points"),
        face = metcor.plot.options$legend$label.fontface
      ),
      #legend.key.size = ggplot2::unit(
      #  metcor.plot.options$legend$size, "lines"
      #),
      legend.margin = ggplot2::margin(
        t = metcor.plot.options$legend$margins$t,
        l = metcor.plot.options$legend$margins$l,
        b = metcor.plot.options$legend$margins$b,
        r = metcor.plot.options$legend$margins$r,
        unit = metcor.plot.options$legend$margins$unit
      ),
      legend.background = ggplot2::element_rect(
        fill = ggplot2::alpha(
          metcor.plot.options$legend$background.color,
          metcor.plot.options$legend$background.alpha
        ),
        size = 0, linetype = 0
      )
    ))

  }


  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                         PREPARING WORLD MAP                         ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################

  # read and prepare the wmap shape for plotting
  if (verbose) {
    message("- preparing world map")
  }
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

  # try to make a grid raster
  # e <- as(raster::extent(-180, 180, -90, 90), "SpatialPolygons") %>%
  #  sf::st_as_sf()


  graticules <- metcor.ne.graticules[["medium"]][["gr_10"]]

  ## trim the raster
  metcor.raster <- terra::trim(metcor.raster)

  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                       WORKING ON CENTER POINT                       ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  if (metcor.plot.options$plot$center.from == "raster") {
    if (verbose) {
      message("- raster is used for mid point calculation, based on coordinate
              distribution peaks")
    }

    # turn raster into tibble
    metcor.raster.tibble <- as_tibble(metcor.raster,
      xy = TRUE,
      na.rm = TRUE
    )

    # find the peak in the lat and lon distribution
    x_center <- stats::density(metcor.raster.tibble$x)$x[which.max(density(metcor.raster.tibble$x)$y)]
    y_center <- stats::density(metcor.raster.tibble$y)$x[which.max(density(metcor.raster.tibble$y)$y)]
    # create a simple feature from this
    zoom_to_xy <- sf::st_transform(sf::st_sfc(sf::st_point(c(x_center, y_center)),
      crs = terra::crs(metcor.raster)
    ),
    crs = 4326
    )
    # get zoom_to coordinates
    zoom_to <- sf::st_coordinates(zoom_to_xy)
  }
  if (metcor.plot.options$plot$center.from == "receptor") {
    if (verbose) {
      message("- receptor is used for mid point calculation")
    }
    if (!identical(receptor, NA)) {
      zoom_to <- c("lon" = receptor$lon[[1]], "lat" = receptor$lat[[1]])
    } else {
      cli::cli_abort(c(
        "{.var receptor} must be provided:",
        "x" = "Need to have receptor details to set the center to the receptor."
      ))
    }
  }
  if (metcor.plot.options$plot$center.from == "manual") {
    if (verbose) {
      message("- manually supplied coordinates are used for mid point")
    }
    zoom_to <- c(
      "lon" = metcor.plot.options$plot$center.point[["lon"]],
      "lat" = metcor.plot.options$plot$center.point[["lat"]]
    )
  }

  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                          HANDLING RECEPTOR                          ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################

  # receptor needs to be a tibble with at least columns named lat and lon with
  # the WGS84 latitude and longitude locations
  if (!identical(receptor, NA)) {
    receptor <- sf::st_transform(sf::st_as_sf(
      x = receptor,
      coords = c("lon", "lat"),
      crs = 4326
    ),
    crs = "ESRI:102017"
    )
  }

  ############################################################################
  ############################################################################
  ###                                                                      ###
  ###                         SETTING ZOOM OPTIONS                         ###
  ###                                                                      ###
  ############################################################################
  ############################################################################

  # Lambert azimuthal equal-area projection around center of interest
  projection_crs <- sprintf(
    "+proj=laea +lon_0=%f +lat_0=%f +x_0=0 +y_0=0
                            +datum=WGS84 +units=m +no_defs",
    zoom_to[1],
    zoom_to[2]
  )

  C <- 40075016.686 # ~ circumference of Earth in meters
  x_span <- C / 2^metcor.plot.options$plot$zoom.level
  y_span <- C / 2^(metcor.plot.options$plot$zoom.level + 1) # also sets aspect ratio

  zoom_to_xy <- sf::st_transform(
    sf::st_sfc(sf::st_point(zoom_to), crs = 4326),
    crs = projection_crs
  )

  disp_window <- sf::st_sfc(
    sf::st_point(sf::st_coordinates(zoom_to_xy - c(x_span / 2, y_span / 2))),
    sf::st_point(sf::st_coordinates(zoom_to_xy + c(x_span / 2, y_span / 2))),
    crs = projection_crs
  )

  if (identical(metcor.plot.options$plot$xlim, NA)) {
    if (verbose) {
      message("   - setting automatic xlim")
    }
    metcor.plot.options$plot$xlim <- sf::st_coordinates(disp_window)[, "X"]
  } else {
    if (verbose) {
      message("   - using provided xlim")
    }
  }

  if (identical(metcor.plot.options$plot$ylim, NA)) {
    if (verbose) {
      message("   - setting automatic ylim")
    }
    metcor.plot.options$plot$ylim <- sf::st_coordinates(disp_window)[, "Y"]
  } else {
    if (verbose) {
      message("   - using provided ylim")
    }
  }

  ###########################################################################
  ###########################################################################
  ###                                                                     ###
  ###                   SETTING RASTER COLORING OPTIONS                   ###
  ###                                                                     ###
  ###########################################################################
  ###########################################################################
  if (metcor.plot.options$raster$type == "discrete") {
    # we need to calculate some things in order to apply discrete colors

    if (identical(metcor.plot.options$raster$discrete.colors, NA)) {
      # apply default coloring
      if (verbose) {
        message("- applying default coloring from ArcMap")
      }
      legend.colors <- c("#2892c7", "#a0c79b", "#fafa64", "#fa8d34", "#e81014")
    }

  # calculate the breaks using classInt
    if (length(metcor.plot.options$raster$discrete.breaks) == 1) {
      if (verbose) {
        message(paste0(
          "- calculating '",
          metcor.plot.options$raster$discrete.breaks,
          "' breaks using non-smoothed raster"
        ))
      }
      metcor.plot.options$raster$discrete.breaks <- classInt::classIntervals(
        var = stats::na.omit(as.vector(values(metcor.raster))),
        n = 5,
        style = metcor.plot.options$raster$discrete.breaks
      )
      metcor.plot.options$raster$discrete.breaks <- metcor.plot.options$raster$discrete.breaks$brks
    }

    # check if breaks match the colors
    if ((length(metcor.plot.options$raster$discrete.breaks) - 1) != length(metcor.plot.options$raster$discrete.colors)) {
      # the provided breaks do not match the colors.
      len.legend.breaks <- length(metcor.plot.options$raster$discrete.breaks) - 1
      len.legend.colors <- length(metcor.plot.options$raster$discrete.colors)
      cli::cli_abort(c(
        "Number of breaks do not match the number of colors:",
        "x" = "There are {len.legend.breaks} breaks and
          {len.legend.colors} colors defined in
          the plot options."
      ))
    }


    # create legend labels, if not provided
    if (identical(metcor.plot.options$legend$labels, NA)) {
      if (verbose) {
        message("- creating legend labels based on the break values")
      }
      for (i in seq_along(metcor.plot.options$raster$discrete.colors)) {
        if (i == 1) {
          metcor.plot.options$legend$labels <- c(
            paste(format(round(
              metcor.plot.options$raster$discrete.breaks[[i]],
              metcor.plot.options$legend$label.digits
            ),
            nsmall = metcor.plot.options$legend$label.digits
            ),
            format(round(
              metcor.plot.options$raster$discrete.breaks[[i + 1]],
              metcor.plot.options$legend$label.digits
            ),
            nsmall = metcor.plot.options$legend$label.digits
            ),
            sep = " - "
            )
          )
        } else {
          metcor.plot.options$legend$labels <- c(
            metcor.plot.options$legend$labels,
            paste(format(round(
              metcor.plot.options$raster$discrete.breaks[[i]],
              metcor.plot.options$legend$label.digits
            ),
            nsmall = metcor.plot.options$legend$label.digits
            ),
            format(round(
              metcor.plot.options$raster$discrete.breaks[[i + 1]],
              metcor.plot.options$legend$label.digits
            ),
            nsmall = metcor.plot.options$legend$label.digits
            ),
            sep = " - "
            )
          )
        }
      }
    }


    # check legend.labels against num colors
    if (length(metcor.plot.options$raster$discrete.colors) != length(metcor.plot.options$legend$labels)) {
      len.legend.labels <- length(metcor.plot.options$legend$labels)
      len.legend.colors <- length(metcor.plot.options$raster$discrete.colors)
      cli::cli_abort(c(
        "Number of labels do not match the number of colors:",
        "x" = "There are {len.legend.labels} labels and
          {len.legend.colors} colors defined in
          the plot options."
      ))
    }
  }

  # apply the smoothing which is for presentation purposes and emulates the
  # ArcMap resampling during display (properties of the raster -> display tab).
  if (metcor.plot.options$raster$smooth.factor > 0) {
    if (verbose) {
      message("- applying display smoothing to raster")
    }
    metcor.raster <- terra::disagg(metcor.raster,
      fact = metcor.plot.options$raster$smooth.factor,
      method = "bilinear"
    )
  }

  if (verbose) {
    message("- preparing raster for plotting")
  }

  # create the default ggplot object
  if (verbose) {
    message("- creating the plot")
  }

  metcor.plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = graticules,
      color = metcor.plot.options$plot$graticules.outline,
      linewidth = metcor.plot.options$plot$graticules.linewidth,
      fill = NA
    ) +
    ggplot2::geom_sf(
      data = world,
      color = metcor.plot.options$plot$world.outline,
      linewidth = metcor.plot.options$plot$world.linewidth,
      fill = metcor.plot.options$plot$world.background
    )


  #################################################################
  ##            Check if we need to plot other layers            ##
  #################################################################
  if (length(metcor.plot.options$raster$layers.before) > 0) {
    # add layers
    if (verbose) {
      message("- adding additional provided layers before raster")
    }
    num.layer <- 0
    for (layer in metcor.plot.options$raster$layers.before) {
      num.layer <- num.layer + 1
      if (verbose) {
        message(paste("    + layer", num.layer))
      } 
      # check that the layer is a sf object
      if ("sf" %in% class(layer$data)) {
        metcor.plot <- metcor.plot +
          ggplot2::geom_sf(
            data = layer$data,
            color = layer$color,
            linewidth = layer$linewidth,
            fill = layer$fill
          )
      } else {
        cli::cli_abort(c(
          "{.var layers.before} must contain Simple Feature (SF) or SpatRaster objects:",
          "x" = "The layers defined in {.var layers.before} can only contain
            objects with class {.cls sf} or {.cls SpatRaster}."
        ))
      }
    }
  }

  if (identical(metcor.plot.options$legend$title, NA)) {
    metcor.plot.options$legend$title <- ""
  }

  # check if we need to apply the gradient
  if (metcor.plot.options$raster$type == "gradient") {
    if (verbose) {
      message("- applying gradient")
    }
    metcor.plot <- metcor.plot +
      tidyterra::geom_spatraster(
        data = metcor.raster,
        ggplot2::aes(fill = lyr.1),
        alpha = metcor.plot.options$raster$alpha
      ) +
      ggplot2::scale_fill_gradientn(
        colours = metcor.plot.options$raster$gradient.color.palette,
        na.value = NA,
        name = metcor.plot.options$legend$title
      )
  } else {
    if (verbose) {
      message("- applying legend breaks and colors")
    }
    ## Discrete rasters
    metcor.raster.discrete <- terra::classify(
      metcor.raster,
      metcor.plot.options$raster$discrete.breaks
    )
    # plot
    metcor.plot <- metcor.plot +
      tidyterra::geom_spatraster(
        data = metcor.raster.discrete,
        ggplot2::aes(fill = lyr.1),
        alpha = metcor.plot.options$raster$alpha
      ) +
      ggplot2::scale_fill_manual(
        values = metcor.plot.options$raster$discrete.colors,
        na.value = NA,
        labels = metcor.plot.options$legend$labels,
        na.translate = FALSE,
        name = metcor.plot.options$legend$title,
        drop = FALSE
      )
  }

  #################################################################
  ##            Check if we need to plot other layers            ##
  #################################################################
  if (length(metcor.plot.options$raster$layers.after) > 0) {
    # add layers
    if (verbose) {
      message("- adding additional provided layers after raster")
    }
    for (layer in metcor.plot.options$raster$layers.after) {
      # check that the layer is a sf object
      if ("sf" %in% class(layer$data)) {
        metcor.plot <- metcor.plot +
          ggplot2::geom_sf(
            data = layer$data,
            color = layer$color,
            linewidth = layer$linewidth,
            fill = layer$fill
          )
      } else {
        cli::cli_abort(c(
          "{.var layers.after} must contain Simple Feature (SF) objects:",
          "x" = "The layers defined in {.var layers.after} can only contain
            objects with class {.cls sf}."
        ))
      }
    }
  }

  # check if we need to add receptors
  if ("sf" %in% class(receptor)) {
    if (verbose) {
      message("- adding receptors")
    }
    metcor.plot <- metcor.plot +
      ggplot2::geom_sf(
        data = receptor,
        shape = metcor.plot.options$receptor$shape,
        size = metcor.plot.options$receptor$size,
        colour = metcor.plot.options$receptor$color,
        fill = metcor.plot.options$receptor$fill,
        stroke = metcor.plot.options$receptor$stroke
      )
  }

  # Apply projection
  metcor.plot <- metcor.plot +
    ggplot2::coord_sf(
      xlim = metcor.plot.options$plot$xlim,
      ylim = metcor.plot.options$plot$ylim,
      expand = FALSE,
      crs = projection_crs,
      datum = projection_crs
    )


  # add scale if needed
  if (metcor.plot.options$plot$show.scale == TRUE) {
    if (verbose) {
      message("- adding scale bar")
    }
    metcor.plot <- metcor.plot +
      ggspatial::annotation_scale(
        location = "br",
        plot_unit = "m",
        text_cex = 1
      )
  }
  
  # add compass if needed
  if (metcor.plot.options$plot$show.compass == TRUE) {
    if (verbose) {
      message("- adding compass")
    }
    if (metcor.plot.options$plot$show.scale == TRUE) {
      # add some padding
      metcor.plot <- metcor.plot +
        ggspatial::annotation_north_arrow(location = "br", 
                                          which_north = metcor.plot.options$plot$compass.which_north, 
                                          style = ggspatial::north_arrow_fancy_orienteering,
                                          pad_x = unit(0.0, "in"), 
                                          pad_y = unit(0.2, "in"))
    } else {
      metcor.plot <- metcor.plot +
        ggspatial::annotation_north_arrow(location = "br", 
                                          which_north = metcor.plot.options$plot$compass.which_north, 
                                          style = ggspatial::north_arrow_fancy_orienteering,
                                          pad_x = unit(0.0, "in"), 
                                          pad_y = unit(0.0, "in"))
    }
  }

  # add annotation if not NA
  if (!identical(metcor.plot.options$annotation$text, NA)) {
    if (verbose) {
      message("- adding annotation")
    }
    # offsets for the annotation. The offset is 2 percent of the total range
    annotate.x <- max(metcor.plot.options$plot$xlim) - ((max(metcor.plot.options$plot$xlim) - min(metcor.plot.options$plot$xlim)) * 0.02)
    annotate.y <- max(metcor.plot.options$plot$ylim) - ((max(metcor.plot.options$plot$ylim) - min(metcor.plot.options$plot$ylim)) * 0.02)

    # add annotation
    metcor.plot <- metcor.plot +
      ggplot2::annotate("text",
        x = annotate.x,
        y = annotate.y,
        label = metcor.plot.options$annotation$text,
        vjust = 1,
        hjust = 1,
        size = metcor.plot.options$annotation$fontsize,
        fontface = metcor.plot.options$annotation$fontface
      )
  }

  if (verbose) {
    message("- applying theme")
  }
  # apply styling
  metcor.plot <- metcor.plot +
    ggplot2::theme_minimal() +
    theme_opts

  if (metcor.plot.options$raster$type == "discrete") {
    metcor.plot <- metcor.plot +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(alpha = 1)))
  }

  # don't show legend if metcor.plot.options$legend$show == FALSE
  if (metcor.plot.options$legend$show == FALSE) {
    if (verbose) {
      message("- removing legend")
    }
    metcor.plot <- metcor.plot +
      list(ggplot2::theme(legend.position = "none"))
  }

  # return plot object
  if (verbose) {
    message("Done!")
  }

  #################
  # output
  #################
  # print(metcor.plot)
  output <- list(
    "plot" = metcor.plot,
    "plot.options" = metcor.plot.options,
    "proj.crs" = projection_crs,
    "origin" = zoom_to,
    "stats" = stats,
    call = match.call()
  )
  class(output) <- "me2tools"
  if (show.plot) {
    plot(metcor.plot)
  }
  invisible(output)
}
