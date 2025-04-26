#' Correct for double counting in m.mass variable.
#'
#' @param F_matrix F matrix
#' @param m.mass name of m.mass variable
#' @param dc_species species to calculate the sum of and subtract from m.mass
#' @param threshold maximum allowed threshold for negative values. Will give
#'   warning if threshold is exceeded.
#' @param neg.rm replace negative corrected m.mass with zero.
#'
#' @return new F matrix with corrected m.mass
#'
#' @noRd
#'
#'
correct_double_counting <- function(F_matrix,
                       m.mass = "m.mass",
                       dc_species,
                       threshold = -1,
                       neg.rm = TRUE) {

  # get only the concentrations as we recalculate the percentages later
  F.matrix <- F_matrix %>% filter(factor_profile == "concentration_of_species")

  # determine which factors we need to modify
  model.runs <- unique(F.matrix[c("model_run")])

  new_F_matrix <- tibble()
  for (run.row in seq(1, nrow(model.runs), 1)) {
    model.run <- model.runs[run.row, ]

    F.subset <- F.matrix %>%
      filter(model_run == model.run$model_run)

    # calculate the species sum
    species.sum <- F.subset %>%
      filter(species %in% dc_species) %>%
      select(-factor_profile) %>%
      select(contains("factor")) %>%
      colSums()

    # subtract species.sum from m.mass
    for (n.factor in seq(1, length(species.sum), 1)) {
      # loop over all factors to subtract the correct amount
      F.subset[F.subset$species == m.mass, ][[paste0("factor_0", n.factor)]] <-
        F.subset[F.subset$species == m.mass, ][[paste0("factor_0", n.factor)]] - as.numeric(species.sum[n.factor])

      if (F.subset[F.subset$species == m.mass, ][[paste0("factor_0", n.factor)]] < 0) {
        if (F.subset[F.subset$species == m.mass, ][[paste0("factor_0", n.factor)]] < threshold) {
          cli::cli_warn(c(
            "Sum of species is larger than {.var {m.mass}}",
            "x" = "A negative value < threshold has been found for factor_0{n.factor}!"
          ))
        }
        # check if neg.rm == TRUE
        if (neg.rm == TRUE) {
          F.subset[F.subset$species == m.mass, ][[paste0("factor_0", n.factor)]] <- 0
        }
      }
    }

    ## create percentage matrixes
    factor.data <- F.subset %>% select(contains("factor"), -factor_profile)
    species.sum <- factor.data
    factor.sum <- factor.data
    for (num.factors in seq(1, ncol(factor.data), 1)) {
      species.sum[, num.factors] <- (factor.data[, num.factors] / rowSums(factor.data)) * 100
      factor.sum[, num.factors] <- (factor.data[, num.factors] / sum(factor.data[, num.factors])) * 100
    }
    ## add additional columns
    species.sum <- dplyr::bind_cols(F.subset %>%
                                      select(model_type,
                                             factor_profile,
                                             model_run,
                                             species) %>%
                                      mutate(factor_profile = "percentage_of_species_sum"),
                                    species.sum)
    factor.sum <- dplyr::bind_cols(F.subset %>%
                                     select(model_type,
                                            factor_profile,
                                            model_run,
                                            species) %>%
                                     mutate(factor_profile = "percentage_of_factor_total"),
                                   factor.sum)
    ## combine to one dataframe
    F.subset <- dplyr::bind_rows(F.subset, species.sum, factor.sum)

    if (nrow(new_F_matrix) == 0) {
      new_F_matrix <- F.subset
    } else {
      new_F_matrix <- dplyr::bind_rows(new_F_matrix, F.subset)
    }
  }

  return(new_F_matrix)
}
