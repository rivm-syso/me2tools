#' Fix Hysplit output for use in MetCor
#'
#' The HYSPLIT generated trajectories output contains the the year in the yyyy
#' format, depending how the HYSPLIT model is set up. The yyyy format is not
#' supported by MetCor, as it is expecting the yy format. In order to allow
#' MetCor to use the calculated trajectories, the year in the data file has to
#' be reset to the yy format. This function is able to change the year format
#' in a file and write the output of HYSPLIT into a different file usable for
#' MetCor.
#'
#' Note: when accidentally supplying correctly formatted files, the output
#' remains unchanged.
#'
#' @param source location of the HYSPLIT output file.
#' @param target location of the MetCor formated Hysplit file.
#'
#' @examples
#' \dontrun{
#' file_list <- list.files("source_folder",
#'                         "hysplit_output.txt",
#'                         full.names = FALSE)
#'
#' for (file in file_list) {
#'   metcor_fix_hysplit(
#'     paste0("source_folder/", file),
#'     paste0("target_folder/", file)
#'   )
#' }
#' }
#'
#' @export
#'
#' @seealso \code{\link{metcor_import}}, \code{\link{metcor_project_raster}}, 
#' \code{\link{metcor_plot}}, \code{\link{metcor_plot_options}}, 
#' \code{\link{metcor_export}}
#'
#'
metcor_fix_hysplit <- function(source, target) {
  hysplit.file <- readr::read_lines(source)
  index_start <- stringr::str_which(hysplit.file, "^\\s*1 BACKWARD OMEGA") + 1
  substr(hysplit.file[index_start], 3, 4) <- "  "
  readr::write_lines(hysplit.file, target)
}
