#' Import MetCor results file into R
#'
#' This function reads the MetCor results file into a spatial raster object
#' (\code{CRS} = \code{EPSG:4326 (WGS 84)}).
#' Currently this is limited to the Northern Hemisphere (x=0:360, y=0:90)
#'
#' @param file MetCor results file
#' @param extent The extent of the grid. Valid values are \dQuote{north} for
#'   the Northern hemisphere, \dQuote{south} for the Southern hemisphere and
#'   \dQuote{worl} for the world.
#' @param type The type of grid, can be \dQuote{metcor} or \dQuote{WGS84}.
#'   With the latter option the grid is transformed to a regular lat/lon grid
#'   within the boundaries of \[-180,180\] for longitude and \[-90,90\] for 
#'   latitude. The \dQuote{metcor} option uses the MetCor transformed grid with
#'   longitude set as \[0,360\].
#' @param na.rm remove columns and rows with only NA values? Default is
#'   \code{TRUE}
#'
#' @return spatial object containing the results from MetCor.
#'
#' @export
#'
#' @seealso \code{\link{metcor_project_raster}}, \code{\link{metcor_plot}},
#' \code{\link{metcor_plot_options}}, \code{\link{metcor_export}},
#' \code{\link{metcor_fix_hysplit}}
#'
#' @importFrom readr read_lines
#' @importFrom terra rast
#' @importFrom terra values
#' @importFrom terra trim
#' @importFrom utils read.table
#' @importFrom cli cli_abort
#' @importFrom dplyr bind_cols
#'
metcor_import <- function(file, 
                          extent = "north",
                          type = "wgs84",
                          na.rm = TRUE) {

  # check if file exists
  if (!file.exists(file)) {
    cli::cli_abort(c(
      "{.var file} must be a valid file:",
      "x" = "File '{file}' does not exists."
    ))
  }
  
  # check if extent = "north", "south", "world"
  if (!tolower(extent) %in% c("north", "south", "world")) {
    cli::cli_abort(
      c(
        "Unkown option for the {.var extent} of the grid:",
        "i" = "The {.var extent} can be 'north', 'south' or 'world'.",
        "x" = "You've provided an unsupported {.var extent}."
      )
    )
  }
  
  # check if type = "metcor", "wgs84"
  if (!tolower(type) %in% c("metcor", "wgs84")) {
    cli::cli_abort(
      c(
        "Unkown option for the {.var type} of the grid:",
        "i" = "The {.var type} can be 'metcor' or 'wgs84'.",
        "x" = "You've provided an unsupported {.var type}."
      )
    )
  }

  # read grid data from file
  metcor.data <- utils::read.table(file,
    header = FALSE,
    sep = "\t",
    skip = 6
  )

  if (sum(is.na(metcor.data[ncol(metcor.data)])) == nrow(metcor.data)) {
    # drop last empty column WHEN it only contains NA (weird writing by MetCor)
    metcor.data <- metcor.data[-ncol(metcor.data)]
  }

  # read the nodata value
  nodata.value <- as.integer(gsub(
    pattern = "nodata_value ",
    replacement = "",
    x = readr::read_lines(file,
      skip = 5,
      n_max = 1
    )
  ))

  # set the values for the nodata value to NA
  metcor.data[metcor.data == nodata.value] <- NA

  # set x and y range based on extent
  x_range <- c(0, 360) # this is always the full range
  if (tolower(extent) == "north") {
    y_range <- c(0, 90)
  }
  if (tolower(extent) == "south") {
    y_range <- c(-90, 0)
  }
  if (tolower(extent) == "world") {
    y_range <- c(-90, 90)
  }
  
  # adjust output if type = wgs84
  if (tolower(type) == "wgs84") {
    metcor.data <- dplyr::bind_cols(metcor.data[,181:360],
                                    metcor.data[,1:180])
    # change ranges for x
    x_range <- x_range - 180
  }
  
  # create an empty SpatRaster with the correct dimensions
  metcor.raster <- terra::rast(
    nrows = nrow(metcor.data),
    ncols = ncol(metcor.data),
    resolution = c(
      (max(x_range)-min(x_range)) / ncol(metcor.data),
      (max(y_range)-min(y_range)) / nrow(metcor.data)
    ), # based on # rows of matrix
    xmin = x_range[1],
    xmax = x_range[2],
    ymin = y_range[1],
    ymax = y_range[2],
    crs = "EPSG:4326"
  )
  gc() # trying to catch those pesky error messages from the garbage collector
  # transform the data to a matrix and add it to the raster
  terra::values(metcor.raster) <- as.matrix(metcor.data)

  # check if the raster has values
  if (identical(terra::global(metcor.raster, min, na.rm = TRUE)[1, 1], NA_real_)) {
    cli::cli_warn(c(
      "The output does not contain any data.",
      "i" = "The cells in the raster are all empty."
    ))
  } else {
    # Trim the raster if needed
    if (na.rm) {
      metcor.raster <- terra::trim(metcor.raster)
    }
  }

  # return SpatRaster
  return(metcor.raster)
}
