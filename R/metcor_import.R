#' Import MetCor results file into R
#'
#' This function reads the MetCor results file into a spatial raster object
#' (\code{CRS} = \code{EPSG:4326 (WGS 84)}).
#' Currently this is limited to the Northern Hemisphere (x=0:360, y=0:90)
#'
#' @param file MetCor results file
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
#'
metcor_import <- function(file, na.rm = TRUE) {
  
  # check if file exists
  if (!file.exists(file)) {
    cli::cli_abort(c(
      "{.var file} must be a valid file:",
      "x" = "File '{file}' does not exists."
    ))
  }

  # read grid data from file
  metcor.data <- utils::read.table(file,
    header = FALSE,
    sep = "\t",
    skip = 6
  )

  # drop last empty column
  metcor.data <- metcor.data[-ncol(metcor.data)]

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

  #######
  # The following code is restricted to the Northern Hemisphere (NH)!
  # TODO: update the code to accept other regions (with a few default options)

  # create an empty SpatRaster with the correct dimensions
  metcor.raster <- terra::rast(
    nrows = nrow(metcor.data),
    ncols = ncol(metcor.data),
    resolution = c(
      360 / ncol(metcor.data),
      90 / nrow(metcor.data)
    ), # based on # rows of NH
    xmin = 0, # based on Northern Hemisphere
    xmax = 360, # based on Northern Hemisphere
    ymin = 0, # based on Northern Hemisphere
    ymax = 90, # based on Northern Hemisphere
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
