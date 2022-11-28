#' Check necessary packages for me2tools
#'
#' Performs a check if all needed packages are actually installed on the system.
#' Will provide messages when missing packages are detected.

#' @return Returns boolean: \code{TRUE} if everything has been found and
#' \code{FALSE} when packages are missing.

#' @export
#'
#' @import stringr
#' @importFrom utils installed.packages
#' @importFrom utils packageDescription
#'
me2tools_check_packages <- function() {
  my_packages <- utils::packageDescription("me2tools")

  my_packages <- stringr::str_split(my_packages$Imports,
                                        pattern = ",\n",
                                        simplify = TRUE)

  not_installed <- my_packages[!(my_packages %in% utils::installed.packages()[ , "Package"])]

  all_installed <- TRUE
  if(length(not_installed) > 0) {
    for (package in not_installed) {
      message(paste0("Package '", package ,"' not installed."))
    }
    all_installed <- FALSE
  }

  return (all_installed)

}


