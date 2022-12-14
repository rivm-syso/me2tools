% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metcor_fix_hysplit.R
\name{metcor_fix_hysplit}
\alias{metcor_fix_hysplit}
\title{Fix Hysplit output for use in MetCor}
\usage{
metcor_fix_hysplit(source, target)
}
\arguments{
\item{source}{location of the HYSPLIT output file.}

\item{target}{location of the MetCor formated Hysplit file.}
}
\description{
The HYSPLIT generated trajectories output contains the the year in the yyyy
format, depending how the HYSPLIT model is set up. The yyyy format is not
supported by MetCor, as it is expecting the yy format. In order to allow
MetCor to use the calculated trajectories, the year in the data file has to
be reset to the yy format. This function is able to change the year format
in a file and write the output of HYSPLIT into a different file usable for
MetCor.
}
\details{
Note: when accidentally supplying correctly formatted files, the output
remains unchanged.
}
\examples{
\dontrun{
file_list <- list.files("source_folder",
                        "hysplit_output.txt",
                        full.names = FALSE)

for (file in file_list) {
  metcor_fix_hysplit(
    paste0("source_folder/", file),
    paste0("target_folder/", file)
  )
}
}

}
\seealso{
\code{\link{metcor_import}}, \code{\link{metcor_project_raster}},
\code{\link{metcor_plot}}, \code{\link{metcor_plot_options}},
\code{\link{metcor_export}}
}
