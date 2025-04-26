#' Check for missing suggested packages
#'
#' This function checks whether the specified packages (typically suggested ones)
#' are installed on the user's system. If any are missing, it displays a helpful
#' message recommending how to install them.
#'
#' @return Invisibly returns a character vector of missing package names.
#'
#' @examples
#' \dontrun{
#' me2_check_suggests()
#' }
#'
#' @export
me2_check_suggests <- function() {

  pkgs <- c(
      "classInt",
      "geomtextpath",
      "ggnewscale",
      "ggpmisc",
      "ggspatial",
      "ggtext",
      "MASS",
      "pals",
      "rnaturalearth",
      "rnaturalearthdata",
      "rnaturalearthhires",
      "sf",
      "terra",
      "tidyterra"
    )

  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing)) {
    message("The following suggested packages are not installed:\n",
            paste(missing, collapse = ", "), "\n\n",
            "Install them with:\n",
            "install.packages(c(", paste(shQuote(missing), collapse = ", "), "))")
  } else {
    message("All suggested packages are installed. ✅")
  }

  invisible(missing)
}