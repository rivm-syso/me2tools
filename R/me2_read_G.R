#' Get the G values from the results file (.txt)
#'
#' Typical use of Multilinear Engine version 2 (ME-2) will provide three
#' different files after a run: \dQuote{.dat}: machine readable format,
#' \dQuote{.rsd}: results for the residuals and \dQuote{.txt}: text file with
#' auxiliary information (i.e. headers). The headers in the text (.txt) file
#' are used in me2tools to split the data into several blocks. In this function
#' the G values are read using user provided header information, denoting the
#' \dQuote{start} and \dQuote{end} lines of the block containing the required
#' data.
#'
#' @param me2_txt_file ME2 output file (.txt), containing the results and
#'  auxiliary information.
#' @param tidy_output Should the output of both G and F be reshaped into tidy
#'   data? Default: FALSE
#' @param block_boundaries A list containing the \dQuote{start} and \dQuote{end}
#'   string used to identify the boundaries of the block containing the values
#'   for the F-matrix.
#' @param dates A vector containing the sample dates for the rows in the
#'   G-matrix. If these dates are outputted in the ME-2 output as the
#'   second column (a column of row numbers being the first), then these values
#'   are used when \code{dates = NA}. The values inside the ME-2 output are
#'   overwritten when a vector if \code{dates} are provided. If this second
#'   column with dates is not available new daily dates, starting from
#'   1970-01-01 will be automatically provided.
#' @param factor_mass This is a vector with the same length as the number of
#'   factor profiles. It contains the total mass in concentration units which is
#'   used to transform the G matrix from unity to concentration units.
#' @param rescale_unity In some cases the averages of the G factors are not
#'   equal to unity. By default a warning is given whenever this is the case.
#'   With this parameter set to \code{TRUE} each factor is multiplied by
#'   1/avg(factor), so that the G factors are scaled to unity again.
#' @param tz Parameter to control the timezone when parameter \dQuote{dates} is
#'   not used. Default: 'Etc/GMT-1'
#'
#' @section Adding dates to G:
#' By far the easiest way to add \dQuote{dates} to G is to provide them as
#' input parameters. The dates can probably be found in the original data
#' input used for ME-2 calculations. For multi time ME-2 application these dates
#' correspond with the lowest time resolution.
#'
#' @section Renaming factor profiles:
#' As the labels of the factor profiles are unknown, this routine reads them as
#' \dQuote{factor_xx}, with xx being a unique number and outputted as a
#' character array. The easiest way to work with this data is by re-coding this
#' vector as \code{factor}. Then the order of the factor profiles and the
#' correct names can be easily set using the following code.
#'
#' ```R
#' mydata$factor <- factor(mydata$factor)
#' mydata$factor <- dplyr::recode_factor(mydata$factor,
#'                                       `factor_01` = "MyFirstName",
#'                                       `factor_02` = "MySecondName",
#'                                       ...
#' )
#' ```
#'
#' Please note that the above will only work when the data is read with the
#' \code{tidy_output = TRUE} setting.
#'
#' @return (tidied) tibble containing G for multiple runs
#'
#' @export
#'
#' @seealso \code{\link{me2_BS_read_G}}, \code{\link{me2_DISP_read_G}},
#' \code{\link{me2_read_all}}, \code{\link{me2_read_dat}}
#'
#' @import cli
#' @import readr
#' @import stringr
#' @import tibble
#' @import dplyr
#'
me2_read_G <- function(me2_txt_file,
                       tidy_output = FALSE,
                       block_boundaries = list("start" = "^\\s*Time factors",
                                               "end" = "^\\s*Composition factors"),
                       dates = NA,
                       factor_mass = NA,
                       rescale_unity = FALSE,
                       tz = "Etc/GMT-1") {

  # check if file exists
  if (!file.exists(me2_txt_file)) {
    cli::cli_abort(c(
      "{.var me2_txt_file} must be a valid file:",
      "x" = "File '{me2_txt_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_txt_file)
  text_filter <- text[text != ""]

  # check if we only have one parameter block
  me2_dupl_block_check(text_filter)

  # sometimes the space for the rownumber is not large enough,
  # leading to ****. In order to read the file we need to replace
  # these with an actual number (-999)
  text_filter <- gsub(" \\*\\*\\*\\*", " -999", text_filter)

  # factor contributions are located between "   Time factors" and
  # "   Composition factors" (note the spaces)
  index_start <- stringr::str_which(text_filter, block_boundaries$start)
  if (length(index_start) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The start header of the data block could not be found in '{me2_txt_file}'.",
      "x" = "Did you provide the correct file?"
    ))
  }
  index_end <- stringr::str_which(text_filter, block_boundaries$end)
  if (length(index_end) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The end header of the data block could not be found in '{me2_txt_file}'.",
      "x" = "Did you provide the correct file?"
    ))
  }

  # the length of the index_start will provide how many runs were performed.
  # we need to get all runs
  g_matrix <- tibble::tibble()
  for (run.number in seq(1, length(index_start), 1)) {
    # get G-values
    g_matrix.tmp <- extract_me2_matrix_between(
      text_filter,
      index_start[run.number],
      index_end[run.number]
    )



    # add dates
    ## Check length of date against length of data!
    if (length(dates) > 1) {
      if (length(dates) == nrow(g_matrix.tmp)) {
        if (!any(grepl("identifier", colnames(g_matrix.tmp)))) {
          g_matrix.tmp <- g_matrix.tmp %>%
            tibble::add_column(
              identifier = dates,
              .before = "factor_01"
            )
        } else {
          g_matrix.tmp <- g_matrix.tmp %>%
            mutate(identifier == dates)
        }
        g_matrix.tmp <- g_matrix.tmp %>%
          mutate(identifier = format(identifier, "%Y-%m-%d %H:%M"))
      } else {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var dates} has a different length ({length(dates)})
           compared to the G-matrix ({nrow(g_matrix.tmp)})",
          "x" = "{.var dates} must have the same length as the number of
          rows in the G-matrix"
        ))
      }
    } else {
      if (!is.na(dates)) {
        cli::cli_abort(c(
          "Different lengths:",
          "i" = "{.var dates} has a different length ({length(dates)})
           compared to the G-matrix ({nrow(g_matrix.tmp)})",
          "x" = "{.var dates} must have the same length as the number of
          rows in the G-matrix or should be set as 'NA'"
        ))
      }
    }


    # get mass_factors
    num_factors <- ncol(g_matrix.tmp %>% select(contains("factor_")))
    if (length(factor_mass) == 1) {
      if (!is.na(factor_mass))  {
        cli::cli_abort(c(
          "Length should be equal to the number of factors:",
          "i" = "{.var factor_mass} has length {length(factor_mass)}.",
          "x" = "Length should be equal to {num_factors}"
        ))
      }
    } else {
      if (length(factor_mass)!=num_factors) {
        cli::cli_abort(c(
          "Length should be equal to the number of factors:",
          "i" = "{.var factor_mass} has length {length(factor_mass)}.",
          "x" = "Length should be equal to {num_factors}"
        ))
      }
    }

    g_matrix.tmp <- g_matrix.tmp %>%
      tidy_me2_contributions(
        factor_mass = factor_mass,
        run_number = run.number,
        rescale_unity = rescale_unity,
        tidy_output = tidy_output,
        tz = tz
      )

    # add run_type
    g_matrix.tmp <- g_matrix.tmp %>%
      tibble::add_column(run_type = "base_run", .after = "model_run")

    if (nrow(g_matrix) > 0) {
      g_matrix <- dplyr::bind_rows(g_matrix, g_matrix.tmp)
    } else {
      g_matrix <- g_matrix.tmp
    }
    rm(g_matrix.tmp)
  }
  return(g_matrix)
}
