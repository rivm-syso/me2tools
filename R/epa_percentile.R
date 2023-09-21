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
    cli::cli_warn(c(
      "{.var x} is not numeric or logical.",
      "x" = "Returning {NA}."
    ))
    return(NA_real_)
  }
  if (!is.numeric(prob) || length(prob) != 1L) {
    cli::cli_abort(c(
      "{.var prob}  must be numeric of length one.",
      "x" = "Please adjust your input."
    ))
  }

  if (prob > 1) {
    cli::cli_warn(c(
      "{.var prob} should have a value between 0 - 1.",
      "x" = "Returning NA."
    ))
    return(NA_real_)
  }
  if (na.rm) {
    x <- x[!is.na(x)]
  }

  n <- length(x)

  L <- ((n + 1) * prob)
  I <- as.integer(L)
  # if we have to few samples (18 or less for P05, I can become 0)
  if (I == 0) {
    I <- 1
    F <- 0
    cli::cli_warn(c(
      "The length of the data (n={n}) is too small to calculate the requested percentile.",
      "x" = "The minimum of the data is used instead."
    ))
  } else {
    F <- L - I
  }

  w1 <- 1 - F
  w2 <- F
  w3 <- 0

  if ((n - I) <= 0) {
    P <- x[n] # grab the maximum value as the best estimate
    cli::cli_warn(c(
      "The length of the data (n={n}) is too small to calculate the requested percentile.",
      "x" = "The maximum of the data is used instead."
    ))
    
  } else if ((n - I) == 1) {
    P <- w1 * x[I] + w2 * x[I + 1]
  } else if ((n - I) >= 2) {
    P <- w1 * x[I] + w2 * x[I + 1] + w3 * x[I + 2]
  }
  return(P)
}
