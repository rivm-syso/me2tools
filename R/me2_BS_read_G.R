#' Get the G values from the BS results file (.txt)
#'
#' Files stored after a BS run are named with a user-specific prefix, shown 
#' here as an asterisk (*). Three output files (*_BS.dat,  *_BS.txt and 
#' *_BS.rsd) are stored after a BS run, and in this function the data provided 
#' in the \dQuote{.txt} file are read. The function \code{me2_BS_read_F}  
#' reads all the factor profiles in the text file and 
#' \code{me2_BS_read_G} (this function) reads all the factor contributions. 
#' Based on the provided \code{corr_threshold} in the \code{me2_BS_read_F} and
#' \code{me2_BS_read_G} functions only factors mapped to base factors with a 
#' correlation larger than \code{corr_threshold} are retained when aggregating 
#' the BS results (i.e., \dQuote{BS_P05}, \dQuote{BS_P95} and 
#' \dQuote{BS_median}). The residuals, stored in the *_BS.rsd file, can be 
#' read using \code{me2_read_residuals}.
#'
#' @param me2_bs_txt_file output file (.txt), containing the results and
#'  auxiliary information for the BS runs.
#' @param corr_threshold This parameter is used to filter out the factors that
#'   have a correlation less to the provided threshold (default = 0.6).
#'   This way only the BS results for factors with a correlation larger or 
#'   equal to the threshold are retained for the aggregated results. The
#'   threshold is also used in the mapping. When the maximum correlation for the
#'   BS factors is below the threshold, the factor is considered to be 
#'   "unmapped".
#' @param tidy_output Should the output be reshaped into tidy data? Default:
#'   FALSE
#' @param block_boundaries A list containing the \dQuote{start} and \dQuote{end}
#'   string used to identify the boundaries of the block containing the values
#'   for the G-matrix. The \dQuote{start_Gmap} and \dQuote{end_Gmap} are used to
#'   identify the boundaries of the block with the G correlations which
#'   are use to produce a mapping table (i.e., which factor has the highest
#'   correlation) that can be used to check swaps between factors.
#' @param dates A vector containing the sample dates for the rows in the
#'   G-matrix. If these dates are outputted in the ME-2 output as the
#'   second column (a column of row numbers being the first), then these values
#'   are used when \code{dates = NA}. The values inside the ME-2 output are
#'   overwritten when a vector of \code{dates} are provided.If this second
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
#'   not used. Default: 'Etc/GMT-1'#'
#'
#' @section Defining the block with data:
#'
#' The \dQuote{.txt} file contains a lot of empty lines, which are stripped
#' immediately after reading the file. This is important, as this ensures that
#' each block of data is contained withing two identifying strings.
#'
#' The function uses these lines to determine the start and end line in which
#' the data lies. As such, the lines between start and end are read and parsed.
#'
#' The identifying strings are provided using the \code{block_boundaries}, which
#' is a list containing the \dQuote{start} and \dQuote{end} string. The string
#' can be preceded by \code{^\\s*} which tells the function that the string can
#' begin with an undefined number of spaces.
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
#' \preformatted{
#' mydata$factor <- factor(mydata$factor)
#' mydata$factor <- dplyr::recode_factor(mydata$factor,
#'                                       `factor_01` = "MyFirstName",
#'                                       `factor_02` = "MySecondName",
#'                                       ...
#' )
#' }
#'
#' Please note that the above will only work when the data is read with the
#' \code{tidy_output = TRUE} setting.
#'
#' @return \code{me2_BS_read_G} returns an object of class ``me2tools''.
#'   The object includes three main components: \code{call}, the command used
#'   to read the data; \code{data}, the G data for each BS run; \code{G_format},
#'   the tidied G data for each BS run where the correlation between the BASE 
#'   and BS factor is larger than the \dQuote{corr_threshold} and the 
#'   \code{G_mapping}, the mapping data based on G correlations of each BS run. 
#'   
#'   The output of the G-values for the BS results is nearly identical to the 
#'   output of the G-values from the base runs. The only difference is that the
#'   BS results have different values in the \dQuote{run_type} column. In this 
#'   case this column contains information of the type of bootstrap.
#'   
#'   The mapping data is constructed with the BASE factors in the columns and
#'   the row factors being the BS factors. Hence, by reading the rows one can
#'   observe swaps from BS factors into different BASE factors. Ideally, a
#'   BS factor should only be mapped to the same BASE factor. When the maximum 
#'   correlation for the BS factors is below the \dQuote{corr_threshold}, the 
#'   factor is considered to be "unmapped".
#'   
#'   If retained, e.g., using \code{output <- me2_BSDISP_read_res(file)}, this 
#'   output can be used to recover the data, reproduce, or undertake further 
#'   analysis.
#'
#'   An me2tools output can be manipulated using a number of generic operations,
#'   including \code{print}, \code{plot} and \code{summary}.
#'
#' @export
#' 
#' @noMd
#'
#' @seealso \code{\link{me2_read_G}}, \code{\link{me2_DISP_read_G}}, 
#' \code{\link{me2_read_all}}, \code{\link{me2_read_dat}}
#'
#'
me2_BS_read_G <- function(me2_bs_txt_file,
                          corr_threshold = 0.6,
                          tidy_output = FALSE,
                          block_boundaries = list(
                            "start" = "^\\s*Factor matrix AA*",
                            "end" = "^\\s*Factor matrix BB",
                            "start_Gmap" = "^\\s*Correlations of G factors:*",
                            "end_Gmap" = "^\\s*Regression matrix T1 of fitted G vs. reference G:*"
                          ),
                          dates = NA,
                          factor_mass = NA,
                          rescale_unity = FALSE,
                          tz = "Etc/GMT-1") {

  # check if file exists
  if (!file.exists(me2_bs_txt_file)) {
    cli::cli_abort(c(
      "{.var me2_bs_txt_file} must be a valid file:",
      "x" = "File '{me2_bs_txt_file}' does not exists."
    ))
  }

  text <- readr::read_lines(me2_bs_txt_file)
  text_filter <- text[text != ""]
  
  # remove '****'
  text_filter <- stringr::str_replace_all(text_filter, '\\*', '9')

  # factor profiles are located between "   Factor matrix AA" and
  # "   Factor matrix BB" (note the spaces)
  index_start <- stringr::str_which(text_filter, block_boundaries$start)
  if (length(index_start) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The start header of the data block could not be found in '{me2_bs_txt_file}'.",
      "x" = "Did you provide the correct file?"
    ))
  }
  index_end <- stringr::str_which(text_filter, block_boundaries$end)
  if (length(index_end) == 0) {
    cli::cli_abort(c(
      "Header not found:",
      "i" = "The end header of the data block could not be found in '{me2_bs_txt_file}'.",
      "x" = "Did you provide the correct file?"
    ))
  }
  
  ## Get the locations of the individual correlations of F for each
  ## BS run in order to determine the mapping of the factor.
  
  index_Gmap_start <- stringr::str_which(text_filter, block_boundaries$start_Gmap)
  index_Gmap_end <- stringr::str_which(text_filter, block_boundaries$end_Gmap)
  # the first block does not have the end string, so we need to add a dummy
  index_Gmap_end <- c(index_Gmap_start[1]+1, index_Gmap_end)
  Gmapping <- tibble()
  
  # add a progress bar
  cli::cli_progress_bar("Reading data", total = length(index_start))

  # the length of the index_start will provide how many runs were performed.
  # we need to get all runs
  g_matrix <- tibble::tibble()
  for (run.number in seq(1, length(index_start), 1)) {
    g_matrix.tmp <- extract_me2_matrix_between(
      text_filter,
      index_start[run.number],
      index_end[run.number],
      headers = FALSE
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
      if (!is.na(factor_mass)) {
        cli::cli_abort(c(
          "Length should be equal to the number of factors:",
          "i" = "{.var factor_mass} has length {length(factor_mass)}.",
          "x" = "Length should be equal to {num_factors}"
        ))
      }
    } else {
      if (length(factor_mass) != num_factors) {
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

    # add column with information
    G.type <- trimws(gsub("Factor matrix AA",
                          "",
                          text_filter[index_start[run.number]]))
    bootstrap.type <- "unknown"
    if (G.type == "Best-fit with original values") {
      bootstrap.type <- "BS_base"
    }
    if (G.type == "Bootstrapped values") {
      bootstrap.type <- "BS"
    }

    if (G.type == "Bootstrapped/pulled values") {
      bootstrap.type <- "BS_pulled"
    }

    g_matrix.tmp <- g_matrix.tmp %>%
      tibble::add_column(run_type = bootstrap.type, .after = "model_run")

    if (nrow(g_matrix) > 0) {
      g_matrix <- dplyr::bind_rows(g_matrix, g_matrix.tmp)
    } else {
      g_matrix <- g_matrix.tmp
    }
    rm(g_matrix.tmp)
    
    ## Create the Gmap output, but always skip the first run
    if (run.number > 1) {
      # read data
      g_corr_run <- extract_me2_matrix_between(
        text,
        index_Gmap_start[run.number],
        index_Gmap_end[run.number],
        headers = FALSE
      )
      # prepare data container
      if (nrow(Gmapping) == 0) {
        Gmapping <- replace(g_corr_run, g_corr_run > -999, 0)
        Gmapping$unmappped <- 0 # add unmapping variable
      }
      for(row_num in seq(from=1, to=nrow(g_corr_run), by=1)){
        Gmap_i <- which.max(g_corr_run[row_num,])
        if (g_corr_run[row_num,Gmap_i] < corr_threshold) {
          # this factor is considered to be unmapped, so update the last column
          Gmapping[row_num, ncol(Gmapping)] <- Gmapping[row_num, ncol(Gmapping)] + 1  
        } else {
          # update a know factor count
          Gmapping[row_num, Gmap_i] <- Gmapping[row_num, Gmap_i] + 1  
        }
      }
    }
    
    # update progress bar
    cli::cli_progress_update()
    
  }
  
  #filter output according to corr_threshold
  factor.correlations <- me2_BS_read_correlations(me2_bs_txt_file = me2_bs_txt_file,
                                                  tidy_output = TRUE)
  
  # Set id variables
  id_variables <- c("model_type", "unit", "model_run", "run_type", "date")
  
  if (tidy_output) {
    g_matrix_threshold <- g_matrix
  } else {
    # Make the table longer, so we can calculate the mid point
    g_matrix_threshold <- g_matrix %>%
      tidyr::pivot_longer(-dplyr::all_of(id_variables), 
                          names_to = "factor",
                          values_to = "value")
  }
  
  g_matrix_threshold <- left_join(g_matrix_threshold, 
                                  factor.correlations %>% 
                                    select(model_run, factor, corr),
                                  by = c("model_run", "factor")) %>% 
    filter(corr >= corr_threshold)
  
  
  output <- list("data" = g_matrix,
                 "G_format" = g_matrix_threshold,
                 "G_mapping" = Gmapping,
                 call = match.call())
  class(output) <- "me2tools"
  
  # close progress bar
  cli::cli_progress_done()
  
  return(output)
}
