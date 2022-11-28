#' Project MetCor raster results
#'
#' This function changes the projection of the imported MetCor ASCII grid from
#' \code{EPSG:4326 (WGS 84)} to a new projection. Currently the only projection
#' supported is the \code{ESRI:102017} projection (North Pole Lambert Azimuthal
#' Equal Area)
#'
#' The projected grid will always have the same resolution as the original grid.
#'
#' @param raster Imported MetCor ASCII grid (see \link{metcor_import})
#' @param projection CRS for the projection (defaults to ESRI:102017)
#' @param method See \code{terra:project} for more details on method.
#' @param na.rm remove columns and rows with only NA values? Default is
#'   \code{TRUE}
#'
#' @return projected raster object
#'
#' @export
#'
#' @seealso \code{\link{metcor_import}}, \code{\link{metcor_plot}}, 
#' \code{\link{metcor_plot_options}}, \code{\link{metcor_export}}, 
#' \code{\link{metcor_fix_hysplit}}
#'
#' @importFrom terra rast
#' @importFrom terra project
#' @importFrom terra trim
#' @importFrom terra res
#'
metcor_project_raster <- function(raster,
                                  projection = "ESRI:102017",
                                  method = "bilinear",
                                  na.rm = TRUE) {
  if (class(raster) != "SpatRaster") {
    cli_abort(c(
      "{.var raster} must be of a class {.cls SpatRaster}:",
      "x" = "You've supplied a {.cls {class(raster)}} raster."
    ))
  }

  if (projection != "ESRI:102017") {
    cli_abort(c(
      "Invalid {.var projection}:",
      "i" = "Only 'ESRI:102017' is currently supported.",
      "x" = "You've supplied: {projection}."
    ))
  }

  # Currently this function only supports the ESRI:102017 projection. Reason
  # for this is because we need to create a SpatRaster to which we need to
  # project on. The raster below is derived from the ArcMap transformation step.

  # get the resolution of the current raster
  raster.res <- terra::res(raster)
  # since it is only for the NH, the ymax = 0 and ymin = -9009965 (= 90 degrees)
  new.x.res <- 9009965 / (90 * raster.res[[1]])
  new.y.res <- 9009965 / (90 * raster.res[[2]])

  proj.raster.source <- terra::rast(
    resolution = c(new.x.res, new.y.res),
    xmin = -9009965,
    xmax = 9009965,
    ymin = -9009965,
    ymax = 9009965,
    crs = projection
  )

  # project old raster onto new raster
  proj.raster <- terra::project(raster, proj.raster.source, method = method)

  if (na.rm) {
    proj.raster <- terra::trim(proj.raster)
  }

  # return result
  return(proj.raster)
}
