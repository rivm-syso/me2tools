# Startup messages
#' @importFrom utils packageVersion
#' @import cli
#' @importFrom crayon green
#' @importFrom crayon red
#' @importFrom crayon blue
#' @importFrom crayon black

.onAttach <- function(lib, pkg, ...){
  
  startup_msg <- paste0("Welcome to ME2tools version ",
                        packageVersion("me2tools"),".\n\n")
  
  
  # check for packages
  ne <- system.file(package = "rnaturalearth")
  nedata <- system.file(package = "rnaturalearthdata")
  nehires <- system.file(package = "rnaturalearthhires")
  
  if(nzchar(ne)) {
    startup_msg <- paste0(startup_msg, 
                          crayon::green(cli::symbol$tick), " ", crayon::blue("rnaturalearth\n"))
  } else {
    startup_msg <- paste0(startup_msg, 
                          crayon::red(cli::symbol$cross), " ", crayon::blue("rnaturalearth\n"))
  }
  
  if(nzchar(nedata)) {
    startup_msg <- paste0(startup_msg, 
                          crayon::green(cli::symbol$tick), " ", crayon::blue("rnaturalearthdata\n"))
  } else {
    startup_msg <- paste0(startup_msg, 
                          crayon::red(cli::symbol$cross), " ", crayon::blue("rnaturalearthdata\n"))
  }
  if(nzchar(nehires)) {
    startup_msg <- paste0(startup_msg, 
                          crayon::green(cli::symbol$tick), " ", crayon::blue("rnaturalearthhires\n"))
  } else {
    startup_msg <- paste0(startup_msg, 
                          crayon::red(cli::symbol$cross), " ", crayon::blue("rnaturalearthhires\n"))
  }
  
  
  startup_msg <- paste0(startup_msg,"\n",
                        crayon::red(cli::symbol$heart), crayon::magenta(" Happy Source Apportioning "), crayon::red(cli::symbol$heart),
                        crayon::black("\n\nTurn this message off using 'suppressPackageStartupMessages(library(me2tools))'"))
  
  packageStartupMessage(startup_msg)
}


################################################################################
## Fix

if (getRversion() >= "2.15.1") {
  
  # What variables are causing issues?
  variables <- c(
    ".", ":=", "!!", "species", "factor_profile", "model_run", "model_type", 
    "identifier", "plot.data", "ymin", "ymax", "y", "x", "xmax", "xmin",
    "middle", "upper", "lower", "n_wllimit", "n_wulimit", "lyr.1", "IDATE",
    "ITIME", "FDATE", "FTIME", "LATR", "LONR", "lat", "lon", "Datestart.",
    "Time", "Dateend.", "Time.1", "Tzone", "X0", "date.start.org", 
    "time.start.org", "date.end.org", "time.end.org", "tzone", "Begin", 
    "Length", "End", "run_type","value", "numeric_x", "DISP_min", "DISP_max",
    "PANEL", "metcor.ne.graticules", "pred", "predict","conf_lwr","conf_upr",
    "CI.linetype", "pred_lwr", "pred_upr", "PI.linetype",  "plot_type", 
    "group_type", ".x", "dummy","run_number_final", "corr"
  )
  
  # Squash the notes
  utils::globalVariables(variables)
  
}