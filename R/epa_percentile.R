#' Percentiles using weighted average approach
#'
#' \code{epa_percentile} calculates the percentiles using the weighted
#' average approach as outlined in the EPA PMF manual.
#'
#' @param x vector for which the percentile needs to be calculated.
#' @param prob probability, ranges from 0 - 1 (i.e. 0.95 calculates the P95)
#' @param na.rm remove NA values
#' @return requested percentile
#'
#' @export
#'
epa_percentile <- function(x, prob, na.rm = TRUE) {
  x <- sort(x)

  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if (!is.numeric(prob) || length(prob) != 1L) {
    stop("'prob' must be numeric of length one")
  }

  if (prob > 1) {
    warning("prob should be between 0 - 1: returning NA")
    return(NA_real_)
  }
  if (na.rm) {
    x <- x[!is.na(x)]
  }

  n <- length(x)

  L <- ((n + 1) * prob)
  I <- as.integer(L)
  F <- L - I

  w1 <- 1 - F
  w2 <- F
  w3 <- 0

  if ((n - I) >= 2) {
    P <- w1 * x[I] + w2 * x[I + 1] + w3 * x[I + 2]
  } else {
    P <- w1 * x[I] + w2 * x[I + 1]
  }
  return(P)
}
