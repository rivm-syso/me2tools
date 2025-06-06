# Startup messages

.onLoad <- function(libname, pkgname) {
  suppressPackageStartupMessages(
    requireNamespace("ggpmisc", quietly = TRUE)
  )
}

.onAttach <- function(lib, pkg, ...){

  core <- data.frame(
    "package" = c(
      "cli",
      "openair",
      "tidyverse",
      "scales"
    )
  )
  core$system_location <- ""
  core$message <- ""
  
  header <- cli::col_black(
    cli::rule(
      left = cli::style_bold("Attaching core packages"),
      right = paste0("me2tools ", packageVersion("me2tools"))
    )
  )
  
  # loop over core packages
  for (i_core in seq_len(nrow(core))) {
    # check packages
    core[i_core,]$system_location <- system.file(package = core[i_core,]$package)
  
    if(nzchar(core[i_core,]$system_location)) {
      core[i_core,]$message <- paste0(cli::col_green(cli::symbol$tick), " ", cli::col_blue(cli::ansi_align(core[i_core,]$package, max(cli::ansi_nchar(core$package)))))
    } else {
      core[i_core,]$message <- paste0(cli::col_red(cli::symbol$cross), " ", cli::col_blue(cli::ansi_align(core[i_core,]$package, max(cli::ansi_nchar(core$package)))))
    }  
  }
  
  if (nrow(core) %% 2 == 1) {
    core <- rbind(core, list("", "", ""))
  }
  col1 <- seq_len(nrow(core) / 2)
  info <- paste0(core[col1,]$message, "     ", core[-col1,]$message)
  
  startup_msg <-paste0(header, "\n", paste(info, collapse = "\n"))
 
  ### suggested
  suggested <- data.frame(
    "package" = c(
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
      "sf",
      "terra",
      "tidyterra"
    )
  )
  suggested$system_location <- ""
  suggested$message <- ""
  
  header <- cli::col_black(
    cli::rule(
      left = cli::style_bold("Suggested packages"),
      right = "")
  )
  
  # loop over suggested packages
  for (i_suggested in seq_len(nrow(suggested))) {
    # check packages
    suggested[i_suggested,]$system_location <- system.file(package = suggested[i_suggested,]$package)
    
    if(nzchar(suggested[i_suggested,]$system_location)) {
      suggested[i_suggested,]$message <- paste0(cli::col_green(cli::symbol$tick), " ", cli::col_blue(cli::ansi_align(suggested[i_suggested,]$package, max(cli::ansi_nchar(suggested$package)))))
    } else {
      suggested[i_suggested,]$message <- paste0(cli::col_red(cli::symbol$cross), " ", cli::col_blue(cli::ansi_align(suggested[i_suggested,]$package, max(cli::ansi_nchar(suggested$package)))))
    }  
  }
  if (nrow(suggested) %% 2 == 1) {
    suggested <- rbind(suggested, list("", "", ""))
  }
  col1 <- seq_len(nrow(suggested) / 2)
  info <- paste0(suggested[col1,]$message, "     ", suggested[-col1,]$message)
  
  suggested_message <- paste0(
    cli::col_cyan(cli::symbol$info), " ",
    cli::col_black("To get the most out of me2tools, we recommend installing any missing packages from the list below.\n  Adding these packages will enable full plotting and analysis functionality of this package.\n  Run `me2_check_suggests()` for more information.")
  )
  
  startup_msg <-paste0(startup_msg, "\n\n", header, "\n", suggested_message, "\n\n", paste(info, collapse = "\n"))
  
  ### suggested
  optional <- data.frame(
    "package" = c(
      "rnaturalearthhires"
    )
  )
  optional$system_location <- ""
  optional$message <- ""
  
  header <- cli::col_black(
    cli::rule(
      left = cli::style_bold("Optional packages"),
      right = "")
  )
  
  # loop over optional packages
  for (i_optional in seq_len(nrow(optional))) {
    # check packages
    optional[i_optional,]$system_location <- system.file(package = optional[i_optional,]$package)
    
    if(nzchar(optional[i_optional,]$system_location)) {
      optional[i_optional,]$message <- paste0(cli::col_green(cli::symbol$tick), " ", cli::col_blue(cli::ansi_align(optional[i_optional,]$package, max(cli::ansi_nchar(optional$package)))))
    } else {
      optional[i_optional,]$message <- paste0(cli::col_red(cli::symbol$cross), " ", cli::col_blue(cli::ansi_align(optional[i_optional,]$package, max(cli::ansi_nchar(optional$package)))))
    }  
  }
  if (nrow(optional) %% 2 == 1) {
    optional <- rbind(optional, list("", "", ""))
  }
  col1 <- seq_len(nrow(optional) / 2)
  info <- paste0(optional[col1,]$message, "     ", optional[-col1,]$message)
  
  optional_message <- paste0(
    cli::col_cyan(cli::symbol$info), " ",
    cli::col_black("Optional packages are packages that can enhance the output by providing more detail (if needed).")
  )
  
  startup_msg <-paste0(startup_msg, "\n\n", header, "\n", optional_message, "\n\n", paste(info, collapse = "\n"))
  
  startup_msg <- paste0(startup_msg,"\n\n",
                        cli::col_red(cli::symbol$heart), 
                        cli::col_magenta(" Happy Source Apportioning "), 
                        cli::col_red(cli::symbol$heart),
                        cli::col_black("\n\nTurn this message off using `suppressPackageStartupMessages(library(me2tools))`."))
  
  
  
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
    "group_type", ".x", "dummy","run_number_final", "corr", "fit", "eq.label",
    "adj.rr.label", "residual_type", "pmf.concentration", "pmf.uncertainty", 
    "sn", "guidance", "Q.Qexp", "dupe", "f_by", "p_by", "eb_by" , "boxplotColor", 
    "group_var", "myColor"
  )
  
  # Squash the notes
  utils::globalVariables(variables)
  
}