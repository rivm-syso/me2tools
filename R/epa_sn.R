#' S/N ratio according to the EPA-PMF V5 manual
#'
#' \code{epa_sn} calculates the S/N ratio using the approach as outlined in 
#' the EPA PMF manual.
#'
#' @param x vector containing the measurements
#' @param x_unc vector containing the measurement uncertainty associated with
#'   the measurements
#' @param na.rm remove NA values
#' @return S/N ratio
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' 
epa_sn <- function(x, x_unc, na.rm = TRUE) {
  
  # note: the s/n is not always the same as calculated by EPA-PMF v5.
  # Currently the cause for this behavior is unknown.
  
  data <- tibble::tibble(x = x,
                 x_unc = x_unc)
  
  if (na.rm) {
    data <- data %>% 
      na.omit()
  }

  data <- data %>% 
    dplyr::mutate(d = dplyr::if_else(x > x_unc, (x-x_unc)/x_unc, 0))
  
  sn <- mean(data$d)

  return(sn)

}
