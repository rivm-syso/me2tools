#' Check if the string is an actual color representation
#' 
#' This function checks if the input, either a string (i.e. "red") or a 
#' character vector contains strings that represent a color. Code is based on
#' information from https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
#' 
#' @param x a string or character vector to test if they represent strings
#'
#' @returns a named logical vector containing TRUE or FALSE
#' @export
#'
#' @examples
is_color <- function(x) {
  vapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  }, logical(1))
}