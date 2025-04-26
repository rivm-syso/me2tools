#' @title Tidy ME-2 factors
#'
#' @description This function tidies the F-matrix to match the PMFR output. It
#' also allows for correction of the m.mass variable to prevent double counting.
#' This is an internal function.
#'
#' @param F_matrix The preliminary output from functions reading F values
#' @param run_number The number associated with the ME-2 run. Default: 1
#' @param tidy_output Should the factor profiles be reshaped into tidy data?
#'   Defaults to \code{FALSE}
#' @param dc_species Sum of the species that should be extracted from missing
#'   mass to prevent double counting.
#'
#' @noRd
#'
#'
tidy_me2_factors <- function(F_matrix,
                             run_number = 1,
                             tidy_output = FALSE,
                             dc_species = NA) {

  # default F has column identifier, but bootstrap has not
  if (!any(grepl("identifier", colnames(F_matrix)))) {
    F_matrix <- F_matrix %>%
      tibble::add_column(
        identifier = paste0(
          "species_",
          sprintf(
            "%02d",
            seq(1, nrow(F_matrix), 1)
          )
        ),
        .before = "factor_01"
      )
  }

  ## add PMFR identifier columns
  F_matrix <- F_matrix %>%
    mutate(identifier = factor(identifier, levels = identifier)) %>%
    dplyr::rename(species = identifier) %>%
    tibble::add_column(model_type = "ME-2", .before = "species") %>%
    tibble::add_column(factor_profile = "concentration_of_species", .before = "species") %>%
    tibble::add_column(model_run = run_number, .before = "species")

  ## create percentage matrixes
  factor.data <- F_matrix %>% select(contains("factor"), -factor_profile)
  species.sum <- factor.data
  factor.sum <- factor.data
  for (num.factors in seq(1, ncol(factor.data), 1)) {
    species.sum[, num.factors] <- (factor.data[, num.factors] / rowSums(factor.data)) * 100
    # remove errors by division by 0
    species.sum[, num.factors][is.nan(as.matrix(species.sum[, num.factors]))] <- 0

    factor.sum[, num.factors] <- (factor.data[, num.factors] / sum(factor.data[, num.factors])) * 100
  }
  ## add additional columns
  species.sum <- dplyr::bind_cols(
    F_matrix %>%
      select(
        model_type,
        factor_profile,
        model_run,
        species
      ) %>%
      mutate(factor_profile = "percentage_of_species_sum"),
    species.sum
  )
  factor.sum <- dplyr::bind_cols(
    F_matrix %>%
      select(
        model_type,
        factor_profile,
        model_run,
        species
      ) %>%
      mutate(factor_profile = "percentage_of_factor_total"),
    factor.sum
  )
  ## combine to one dataframe
  F_matrix <- dplyr::bind_rows(F_matrix, species.sum, factor.sum)

  if (!is.na(dc_species[1])) {
    # remove double counting!
    F_matrix <- correct_double_counting(F_matrix,
      m.mass = "m.mass",
      dc_species,
      threshold = -1,
      neg.rm = TRUE
    )
  }

  if(tidy_output) {
    # Set id variables
    id_variables <- c("factor_profile", "model_run", "species", "model_type")

    if ("run_type" %in% names(F_matrix)) {
      id_variables <- c(id_variables, "run_type")
    }

    # Make the table longer
    F_matrix <- F_matrix %>%
      tidyr::pivot_longer(-dplyr::all_of(id_variables), names_to = "factor") %>%
      arrange(factor,
              factor_profile,
              species)

  }

  return(F_matrix)
}
