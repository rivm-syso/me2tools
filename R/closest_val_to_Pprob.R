#' Translate a value to the closest percentile probability in a vector
#'
#' The polarPlot function from openair requires a probability between 0 and 100.
#' When multiple sites are evaluated it is best to keep the cut-off value the
#' same. However, contributions for sites are different, therefore the
#' percentiles also differ. To keep this value the same, this function can be
#' used to find the corresponding probability in other vectors.
#'
#' @param x vector with values
#' @param value cut off value for which the closest percentile in vector x needs
#' to be determined
#'
#' @return probability [0,1] which can be used in polarPlot (multiply by 100)
#'
#' @export
#'
closest_val_to_Pprob <- function(x, value) {
  return ((which.min(abs(sort(x)-value)))/length(x))
}
