#' Match factors between base and source based on Pearson correlation
#'
#' With this function the Pearson correlation coefficient between the source
#' and base F or G matrix can be calculated. It can be used to assign a
#' bootstrap (BS) run factor to a base factor and observe the number of
#' unmapped factors. Or it can be used to map the results of any other run
#' towards the selected base case with the lowest Q or Qmain contribution.
#'
#' @param base base F.matrix or base G.matrix. This is only one matrix from
#'   a single run.
#' @param source source F.matrix or source G.matrix. This function can only
#'   handle comparison between one other matrix from a single run.
#' @param F.profile where should the correlation be calculated on when using
#'   the F matrix? Defaults to \code{percentage_of_species_sum}
#' @param G.unit where should the correlation be calculated on when using the
#'   G matrix? Defaults to \code{normalised}
#' @param corr.threshold value for the Pearson correlation threshold, defaults
#'   to 0.6 as stipulated by the EPA-PMF 5.0 manual.
#' @param extended.output provide extended output (besides the factor output
#'   the correlation matrix is also provided when this value is set to
#'   \code{TRUE})
#'
#' @return factor output (named vector). If \code{extended.output} is
#'   \code{TRUE} a list is returned with the factor output and the correlation
#'   matrix.
#'
#' Note, in order to use this function to assess the mapping of multiple sources
#' for F or G matrix, the user should use this function inside a loop.
#'
#' The function can also be used to assess the correlations of F and G factors
#' with partly overlapping species (F) or samples (G). If different amount of
#' species or samples are found, these are added to the matrices who do not have
#' them and their contribution is set to 0. If the user wishes to compare only
#' those species or samples present in both base and source, the input should
#' be modified as such.
#'
#' @export
#'
match_factor_order <- function(base,
                               source,
                               F.profile = "percentage_of_species_sum",
                               G.unit = "normalised",
                               corr.threshold = 0.6,
                               extended.output = FALSE) {

  # base is the leading order, the source is matched to the base

  # check if it is F or G that is presented.

  if ("factor_profile" %in% names(base)) {
    presented.matrix = "F.matrix"
  } else if("unit" %in% names(base)) {
    presented.matrix = "G.matrix"
  } else {
    cli::cli_abort(c(
      "{.var base} is unknown:",
      "i" = "Unable to determine if {.var base} is the F or G matrix.",
      "x" = "Did you provide the full tibble?"
    ))
  }

  # First work with F
  if (presented.matrix == "F.matrix") {
    base <- base  %>%
      filter(factor_profile == !!F.profile) %>%
      select_if(~sum(!is.na(.)) > 0)

    source <- source  %>%
      filter(factor_profile == !!F.profile) %>%
      select_if(~sum(!is.na(.)) > 0)

    # check if there are difference in species.
    diff.base.source <- setdiff(base$species, source$species)
    if(length(diff.base.source > 0)) {
      # there are species in base.F who are not in source.F
      for(item.num in seq(1, length(diff.base.source),1)) {
        # grab the first row from source.F
        new_row <- source[1,]
        new_row$species <- diff.base.source[item.num]
        for (factor.col in seq(6,ncol(source),1)) {
          new_row[factor.col] <- 0  # set all contributions to zero
        }
        source <- dplyr::add_row(source, new_row)
        message(paste0("Added species ", diff.base.source[item.num], " to the source matrix with contributions = 0"))
      }
    }

    diff.source.base <- setdiff(source$species, base$species)
    if(length(diff.source.base > 0)) {
      # there are species in source.F who are not in base.F
      for(item.num in seq(1, length(diff.source.base),1)) {
        # grab the first row from source.F
        new_row <- base[1,]
        new_row$species <- diff.source.base[item.num]
        for (factor.col in seq(6,ncol(base),1)) {
          new_row[factor.col] <- 0  # set all contributions to zero
        }
        base <- dplyr::add_row(base, new_row)
        message(paste0("Added species ", diff.base.source[item.num], " to the base matrix with contributions = 0"))
      }
    }

    # arrange by species
    base <- base %>%
      arrange(species)
    source <- source %>%
      arrange(species)

  }


  if (presented.matrix == "G.matrix") {
    base <- base  %>%
      filter(unit == !!G.unit) %>%
      select_if(~sum(!is.na(.)) > 0)

    source <- source  %>%
      filter(unit == !!G.unit) %>%
      select_if(~sum(!is.na(.)) > 0)

    # check if there are difference in species.
    diff.base.source <- setdiff(base$date, source$date)
    if(length(diff.base.source > 0)) {
      # there are dates in base.F who are not in source.F
      for(item.num in seq(1, length(diff.base.source),1)) {
        # grab the first row from source.F
        new_row <- source[1,]
        new_row$date <- diff.base.source[item.num]
        for (factor.col in seq(6,ncol(source),1)) {
          new_row[factor.col] <- 0  # set all contributions to zero
        }
        source <- dplyr::add_row(source, new_row)
        message(paste0("Added date ", diff.base.source[item.num], " to the source matrix with contributions = 0"))
      }
    }

    diff.source.base <- setdiff(source$date, base$date)
    if(length(diff.source.base > 0)) {
      # there are dates in source who are not in base
      for(item.num in seq(1, length(diff.source.base),1)) {
        # grab the first row from source.F
        new_row <- base[1,]
        new_row$date <- diff.source.base[item.num]
        for (factor.col in seq(6,ncol(base),1)) {
          new_row[factor.col] <- 0  # set all contributions to zero
        }
        base <- dplyr::add_row(base, new_row)
        message(paste0("Added date ", diff.base.source[item.num], " to the base matrix with contributions = 0"))
      }
    }

    # arrange by date
    base <- base %>% arrange(date)
    source <- source %>% arrange(date)
  }

  # rename the source factors
  names(source)[6:ncol(source)] <- paste0(names(source)[6:ncol(source)], "_source")


  # calculate the correlations
  cor_matrix <- cor(source %>%
                      select(-contains("factor_profile")) %>%
                      select(contains("factor_")) %>%
                      select_if(~sum(!is.na(.)) > 0),
                    base %>%
                      select(-contains("factor_profile")) %>%
                      select(contains("factor_")) %>%
                      select_if(~sum(!is.na(.)) > 0)
  )

  # create output
  factor.output <- base[1,6:ncol(base)]
  for (factor.name in names(factor.output)) {
    factor.output[[factor.name]] <- paste0("unmapped_", str_sub(factor.name,-2,-1)) #default setting
  }

  # loop over the cor_matrix, check if the correlation > threshold and store this is the case.
  for (n.col in seq(1,ncol(cor_matrix), 1)) {
    highest.cor <- -999;
    highest.cor.row <- 0;
    for (n.row in seq(1,nrow(cor_matrix), 1)) {
      # first check against threshold
      if (cor_matrix[n.row,n.col] > corr.threshold) {
        # check if it is higher than the highest found
        if (cor_matrix[n.row,n.col] >= highest.cor) {
          highest.cor <- cor_matrix[n.row,n.col]
          highest.cor.row <- n.row
        }
      }
    }
    if (highest.cor.row > 0) {
      factor.output[[names(factor.output)[n.col]]]<- gsub("_source", "", rownames(cor_matrix)[highest.cor.row])
    }
  }
  if (extended.output) {
    return(list("factor.output" = factor.output, "cor.matrix" = cor_matrix))
  } else {
    return(factor.output)
  }
}
