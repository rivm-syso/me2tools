#' Add existing species
#' 
#' This function add species to the F-matrix when they are present in
#' the output or when the user provides them. If they are present and
#' overwritten with the user provided species, flags are set that can
#' be used to inform the user
#'
#' @param f_matrix.tmp The temporary tibble holding the F_matrix values
#' @param species A vector containing the names of the species for the rows in
#'   the F-matrix. If these species name are outputted in the ME-2 output as the
#'   second column (a column of row numbers being the first), then these values
#'   are used when \code{species = NA}. If this second column with names is not
#'   available all species are named as \dQuote{species_xx}, with xx being an
#'   unique number starting at 1.
#'   
#' @return List object with the updated \dQuote{f_matrix.tmp} and 
#'   \dQuote{flags}. The \dQuote{flags} is a list containing two items to
#'   denote if species were already present and if they were overwritten by the
#'   user specified species.
#'   
#' @noRd
#'

add_existing_species <- function(f_matrix.tmp,
                                  species = NA) {

  # check if the first column happens to be a character
  # The user might not be aware that they have species in there
  species.present <- FALSE
  species.overwrite <- FALSE
  if ("character" %in% class(f_matrix.tmp[[1]])) {
    # fix it and warn the user
    names(f_matrix.tmp)[2:ncol(f_matrix.tmp)] <- names(f_matrix.tmp)[1:ncol(f_matrix.tmp)-1]
    names(f_matrix.tmp)[1] <- "identifier"
    species.present <- TRUE
  }
  
  # check if we need to add species
  ## Check length of species against length of data!
  if (length(species) > 1) {
    if (length(species) == nrow(f_matrix.tmp)) {
      if (!any(grepl("identifier", colnames(f_matrix.tmp)))) {
        f_matrix.tmp <- f_matrix.tmp %>%
          tibble::add_column(
            identifier = species,
            .before = "factor_01"
          )
      } else {
        f_matrix.tmp <- f_matrix.tmp %>%
          dplyr::mutate(identifier == species)
        if (species.present) {
          species.overwrite <- TRUE
        }
      }
    } else {
      cli::cli_abort(c(
        "Different lengths:",
        "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(f_matrix.tmp)})",
        "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix"
      ))
    }
  } else {
    if (!is.na(species)) {
      cli::cli_abort(c(
        "Different lengths:",
        "i" = "{.var species} has a different length ({length(species)})
           compared to the F-matrix ({nrow(f_matrix.tmp)})",
        "x" = "{.var species} must have the same length as the number of
          rows in the F-matrix or should be set to 'NA'"
      ))
    }
  }
  
  output <- list("f_matrix.tmp" = f_matrix.tmp,
                 "flags" = list("species.present" = species.present,
                                "species.overwrite" = species.overwrite))
  
  return(output)
}